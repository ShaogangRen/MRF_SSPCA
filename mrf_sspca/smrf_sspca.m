function [D,S,W,A,status] = smrf_sspca(X,Y,D,S,W,A,param,hMRF)
%% function for mrf structured sparse PCA

updateS = 1;
cntST = 0;
status = 0;

%%%%%%%%%%%%%%%%%%  step1: get the structure  %%%%%%%%%%%%%%%%%%%%
for iter = 1: param.MaxIte
  oldD = D; 
  oldS = S;
    
 %%%%% compute A %%%%%%%%
     [sD, state]= fix_SD(S,D);
      if state == 1
          status = 1;
          fprintf('Trivial sD was generated, stop! \n');
          break;
      end
     
      sDW = [sD;W]; XY = [X;Y];
    % A = inv(sD'*sD + Lambda6)*(X'*sD)';
      for aiter =  1 : param.AMaxIte
             oldA = A;
             permK = randperm(param.K);
             for tk = 1:param.K
                 k = permK(tk);
                 Aa = (sDW(:,k))'*(XY - (sDW*A - sDW(:,k)*A(k,:)))/((sDW(:,k))'*sDW(:,k));
                 A(k,:) = Aa/sqrt(Aa*Aa');              
             end

            diff_A = norm(A - oldA, 'fro')/norm(oldA, 'fro');
            if diff_A < param.AOpt
                break;
            end 
      end
      
      
%%%%%%%%% compute D S %%%%%%
    for dsIte = 1: param.DSMaxIte
      oldDs = D;
    %%% update D 
      tX = X + ((1-S).*D)*A;
      D = Subgradient_D_V_Eight(tX,A,param);
    %%% update S
      [S,energy] = mrfSegment(X,S,D,A,hMRF,param);  

       diff_Ds = norm(D - oldDs, 'fro')/norm(oldDs,'fro');
       if diff_Ds < param.DsOpt
           break;
       end 
    end
    
%%%%%%% Compute W %%%%%
     W = Subgradient_W_V_Eight(Y,A,param);
%%%   
    if  sum(sum(abs(D))) == 0
        fprintf('\n D all zeros! \n');
        status = 1;
        break;
    end

   diff_D = norm(D - oldD, 'fro')/norm(oldD,'fro');
   diff_S = norm(S - oldS, 'fro')/norm(oldS,'fro');
   sD = S.*D;
   XErr = norm(X - sD*A,'fro')/norm(X,'fro');
    
   
    if param.DebugPrint
            fprintf('%s \n', ...
                        [ sprintf('Iter = \t'),int2str(iter), ...
                          sprintf('\t diff_D = \t'),  num2str(diff_D),...
                          sprintf('\t diff_S = \t'),  num2str(diff_S), ...
                          sprintf('\t XErr = \t'),  num2str(XErr)...
                          ] ); 
    end

    if isnan(diff_D)
        fprintf('\n NaN error! \n');
        status = 2;
        break;
    end
    
    if diff_D < param.DOpt  
        break;
    end 
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sD, state]= fix_SD(S,D)
sD = S.*D;
K = size(D,2);
state = 0;
cnt = 0;
%%%%%%%  sD(:,k) should not be all zero
for k = 1:K
     if sum(abs(sD(:,k))) == 0
        cnt = cnt + 1;
        if sum(abs(D(:,k))) ~= 0
            fprintf('Warning! sDk is replaced with Dk!\n');
            sD(:,k) = D(:,k);
        elseif sum(abs(S(:,k)))~= 0
            % use the sum of rest dictionaries
            fprintf('Warning! sDk is replaced with ssDk!\n');
            ssD =  sum(sD,1)/(K-1);
            sD(S(:,k),k) = ssD(S(:,k));
        else
            fprintf('Warning! sDk is replaced with random values!\n');
            sD(:,k) = 2*(rand(1,size(D,1)) - 0.5);
        end
     end
end
     
if cnt > K*0.5
  state = 1;
end

    


