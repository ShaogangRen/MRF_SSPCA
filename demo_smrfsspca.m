
clear all;
close all;

%addpath('..\KSVD_Matlab_ToolBox');
addpath('gco-v3.0/matlab');
addpath('gco-v3.0/matlab/bin');
addpath('smrf_sspca/');

dataFile = 'data/ARtrainS0.2';
load([ dataFile '.mat']);


%==== model parameters ===
param.lam1 = 40;   % sparsity control lamda1, S_L1
param.lam2 = 30;    % smooth control, |s_ab - s_cd|
param.lam3 = 700;   %% lambda for  ||D||_1
param.lam4 = 10;   %% lambda for  ||W||


param.K = 36;
param.imheight = height;
param.imwidth = width;
param.impixelN = param.imheight*param.imwidth;
param.sampleN = size(X,2);

%%%%%%%%%%%%%%%%% optimization %%%%%%%%%%%%%%%%%%%%%%%

param.MaxIte = 60;
param.DSMaxIte = 80; %Dk Sk combination
param.DsubMaxIte =  200;  %D subgradient
param.SMaxIte = 3;   % graph cut
param.AMaxIte  = 5;
param.rMaxIte = 100;
param.rAMaxIte  = 30;
param.rDMaxIte = 100;


param.AOpt = 0.01; % inner A
param.DsubOpt = 0.01; % subgradient D
param.SkOpt = 0.001; % Sk Dk  combination 
param.DsOpt = 0.1; % D S
param.DOpt = 0.06; % outer D
%%%%% Revise
param.rDOpt = 0.01;
param.rAOpt = 0.01;
param.rdOpt = 0.01;


%% debug control
param.DebugPrint = 1;

fileNML = ['Dict_K' num2str(param.K) 'LOne.mat'];

%%%%%%%%%%%%% Algorithm start %%%%%%%%%%%%%%%%%%
[iD,iS,iW,iA] = InitialDSWA(X,Y,param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Initialize Graph...');
hMRF = GCO_Create(numel(iS(:,1)),2);
GCO_SetSmoothCost( hMRF, [0 1;1 0] );
AdjMatrix = getAdj([param.imheight,param.imwidth]);
GCO_SetNeighbors( hMRF, param.lam2 * AdjMatrix );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,S,A] = smrf_sspca(X,Y,iD,iS,iW,iA,param,hMRF);
GCO_Delete(hMRF);
%%%%%%%%%%%%%%%%%%%%%

fprintf(['File saved to ' fileNML '\n']);
save(fileNML1,'X','D','S','A');



