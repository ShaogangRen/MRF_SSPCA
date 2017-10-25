function [D,S,A] = mrf_sspca(X,D,S,A,param,hMRF)
%% function for mrf structured sparse PCA

updateS = 1;
cntST = 0;

for iter = 1: param.MaxIte
  oldD = D; 
  oldS = S;
    
    %%%%% compute A %%%%%%%%
     sD = S.*D;

      for akiter =  1 : param.AMaxIte
             oldA = A;
             permK = randperm(param.K);
             for tk = 1:param.K
                 k = permK(tk);
            
               %%%%%%%  sD(:,k) should not be all zero
                 if sum(abs(sD(:,k))) == 0
                    if sum(abs(D(:,k))) ~= 0
                      fprintf('Warning! sDk is replaced with Dk!\n');
                        sD(:,k) = D(:,k);
                    elseif sum(abs(S(:,k)))~= 0
                        fprintf('Warning! sDk is replaced with ssDk!\n');
                        ssD =  sum(sD,1);
                        sD(S(:,k),k) = ssD(S(:,k));
                    else
                      fprintf('Warning! sDk is replaced with random values!\n');
                      sD(:,k) = 2*(rand(1,size(D,1)) - 0.5);
                    end
                 end
                %%%%%%% end  %%%%
                
                 Aa = (sD(:,k))'*(X - (sD*A - sD(:,k)*A(k,:)))/((sD(:,k))'*sD(:,k));
                 A(k,:) = Aa/sqrt(Aa*Aa');              
             end

            diff_A = norm(A - oldA, 'fro')/norm(oldA, 'fro');
            if diff_A < param.AOpt
                break;
            end 
      end
  
      
%% compute D S %%%%%% 
    for dsIte = 1: param.DSMaxIte
      oldDs = D;
    %%% update D 
      tX = X + ((1-S).*D)*A;
      D = Subgradient_D(tX,A,param);
    %%% update S
      [S,energy] = mrfSegment(X,S,D,A,hMRF,param);  

       diff_Ds = norm(D - oldDs, 'fro')/norm(oldDs,'fro');
       if diff_Ds < param.DsOpt
           break;
       end 
    end
%%%%%%%%%%%%%%%%%%%%%%

   diff_D = norm(D - oldD, 'fro')/norm(oldD,'fro');
   diff_S = norm(S - oldS, 'fro')/norm(oldS,'fro');
   
    if param.debug_out 
            fprintf('%s \n', ...
                        [ sprintf('Iter = \t'),int2str(iter), ...
                          sprintf('\t diff_D = \t'),  num2str(diff_D),...
                          sprintf('\t diff_S = \t'),  num2str(diff_S), ...
                          ] ); 
    end
    
    
    if diff_D < param.DOpt  || isnan(diff_D)
        break;
    end 
end
