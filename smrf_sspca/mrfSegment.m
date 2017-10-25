function [S, en] = mrfSegment(X,S,D,A,hMRF,mrfpm)

%MaxCost = 10000000

npx = size(X,1);% number of image pixels
N = size(X,2);% number of samples


%%#######################################################
D2 = D.*D;
A2 = A.*A;
XAt = X*A';
DXAt = D.*XAt;

A2sum = sum(A2, 2);

preS = 0;

for itr = 1:mrfpm.SMaxIte
   
    
	en= [0.0 0.0 0.0];
    for k = 1 : mrfpm.K
        
    mrf_k = hMRF;
    
   % coeficient for s = 0
    s0c = 2*DXAt(:,k);  
        
   % compute coeficient for s=1
    sD = S.*D;
    sDA = sD*A;

%=================================   
    sDAk = sDA - sD(:,k)*A(k,:);
    sDAkAt = sDAk*A';
    s1c = D(:,k).*(2*sDAkAt(:,k)) + D2(:,k)*A2sum(k);
%%============================      
    s1c = s1c + mrfpm.lam1;   %%%!!!
    
    GCO_SetDataCost(mrf_k, [s0c, s1c]' );
    GCO_Expansion(mrf_k);

    [ek,dk,sk] = GCO_ComputeEnergy(mrf_k);
    en(1) = en(1) + ek;
    en(2) = en(2) + dk;
    en(3) = en(3) + sk;

    %en= en + 1.0*[ek dk sk];
    S(:,k) = reshape(GCO_GetLabeling(mrf_k)==2,size(S(:,k)));
    end

    S_err = norm(S - preS,'fro');

    if(S_err <= 2)
        break;
    end
    
    preS = S;
end

%fprintf('\n');
%%#############################################################

function W = getAdj(dsize)
numSites = prod(dsize);
id1 = [1:numSites, 1:numSites];
%id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 1+1:numSites+1,...
        1+dsize(1):numSites+dsize(1)]; %...
        %1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
value = ones(1,2*numSites);
W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);







