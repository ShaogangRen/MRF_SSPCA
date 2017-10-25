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

