function D = Subgradient_D(X,A,param)
%shooting method to minimize the following formula
%min_{D} ||x-DA||^2_F+\lambda2||D||

G = -2*X*A';
Aa = A*A';
%Qq = -param.lam2*(Q+Q');

dlam = 1;
DLam = dlam* eye(size(A,1));

D = (inv((A*A'+DLam))*(X*A')')';

iter = 0;
[r,c] = size(D);
converged = 0;

%%%%%%%%%%%%%%%%%%%%%%%
while ~converged && (iter < param.DsubMaxIte)
 subIte = iter; 
 oldDsub = D;
  
  for i = 1:r
      for j = 1:c
        alij = 2*Aa(j,j) ; %+ Qq(i,i);  
        betij = -(2*D(i,:)*Aa(:,j) - 2*D(i,j)*Aa(j,j) + G(i,j));
      %  Qq(i,:)*D(:,j)-Qq(i,i)*D(i,j) +
      %  betij = 0;
      
        if betij < - param.lam1
         D(i,j) =  (betij+ param.lam1)/ alij; 
        elseif betij > param.lam1
         D(i,j) =  (betij - param.lam1)/ alij;    
        else
         D(i,j) = 0; 
        end
      end
  end
  
  diff_Dsub = norm(D - oldDsub, 'fro')/norm(oldDsub,'fro');
  converged = ( diff_Dsub < param.DsubOpt);
  iter = iter + 1;
end

















