function [D,S,W,A] = InitialDSWA(X,Y,param)
% Initialize D, A and S.
%
D = rand(param.impixelN,param.K);
D  = 2*(D - 0.5);
%A = rand(param.K, param.sampleN);

A = rand(param.K,size(X,2));

A = 2*(A - 0.5);

%%%%
W = rand(size(Y,1),param.K);
W = 2*(W - 0.5);

%%% randomly generate S %%%%%%%
S = ones(size(X,1),param.K);
% S = [];
%  for ik =1:param.K
%      tmpidx = randperm(param.impixelN);
%      tmpv = zeros(param.impixelN,1);
%      tmpv(tmpidx(1:ceil(param.impixelN*param.percent))) = 1;
%      S = [S,tmpv];
%  end







