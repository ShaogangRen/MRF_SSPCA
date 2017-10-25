
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRF Structured sparse PCA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%addpath('..\KSVD_Matlab_ToolBox');
addpath('gco-v3.0\matlab');
addpath('mrf_sspca\');
%addpath(genpath('..\L1NormCode\SLEP_package_4.1\SLEP'));

dataFile = 'data\sim';
load([ dataFile '.mat']);

%%%%%%%%%%%%%%%%% Parameter Setting %%%%%%%%%%%%%%%%%%%%%%%
param.K = 3;  %% component number
param.imheight = height;
param.imwidth = width;
param.impixelN = param.imheight*param.imwidth;
param.sampleN = size(X,2);


%====  mrf parameters ===
param.lam2 = 0.3;    % smooth control, |s_ab - s_cd|
param.lam1 = 40;    % sparsity control lamda1, S_L1
param.lam3 = 300    % lambda for  ||D||_1

%====== optimization ====
param.MaxIte = 190;
param.DSMaxIte = 180; %Dk Sk combination
param.DsubMaxIte =  200;  %D subgradient
param.SMaxIte = 3;   % graph cut
param.AMaxIte  = 5;


param.AOpt = 0.01; % inner A
param.DsubOpt = 0.01; % subgradient D
param.SkOpt = 0.001; % Sk Dk  combination 
param.DsOpt = 0.1; % D S
param.DOpt = 0.027; % outer D

%% debug control
param.debugDA = 0;
param.debugS  = 0;
param.debug_out = 1;


[iD,iA,iS] = InitialDAS(X,param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Initialize Graph...');
hMRF = GCO_Create(numel(iS(:,1)),2);
GCO_SetSmoothCost( hMRF, [0 1;1 0] );
AdjMatrix = getAdj([param.imheight,param.imwidth]);
GCO_SetNeighbors( hMRF, param.lam2 * AdjMatrix );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D,S,A] = mrf_sspca(X,iD,iS,iA,param,hMRF);
GCO_Delete(hMRF);
 

fileNM = ['Dict_K' num2str(param.K) '.mat'];
fprintf(['File saved to ' fileNM '\n']);
save(fileNM,'X','D','S','A');




DS = S.*D;

figure('Name','DS1');
imagesc( reshape( DS(:,1), height, width ) )
pause

figure('Name','DS2');
imagesc( reshape( DS(:,2), height, width ) )
pause

figure('Name','DS3');
imagesc( reshape( DS(:,3), height, width ) )
pause



