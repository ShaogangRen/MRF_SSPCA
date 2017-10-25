function [D,A,S] = InitialDAS(X,param)
% Initialize D, A and S.
%
D = rand(param.impixelN,param.K);
D  = 2*(D - 0.5);
%A = rand(param.K, param.sampleN);

A = rand(param.K,size(X,2));

A = 2*(A - 0.5);
  
%%% randomly generate S %%%%%%%
S = ones(size(X,1),param.K);







