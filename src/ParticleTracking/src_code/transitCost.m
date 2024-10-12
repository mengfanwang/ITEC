function [neibIdx, Cij] = transitCost(xCoord,yCoord,zCoord,g,q)
% for transition cost
neibIdx = cell(g.particleNum, 1);
%neiFm = cell(length(movieInfo.xCoord), 1);
Cij = cell(g.particleNum, 1);
curLen = 0;
for i=1:q.timepts
    disp(i);
    for j=i+1:i+g.k
        if j>q.timepts
            break;
        end
        D = createDistanceMatrix([xCoord{i}, yCoord{i}, zCoord{i}],[xCoord{j}, yCoord{j}, zCoord{j}]);
        bsIdxLen = length(cat(1, xCoord{1:j-1}));
        for nn = 1:size(D,1)
            tmpIdx = find(D(nn,:) < g.maxDist);
            if ~isempty(tmpIdx)
                neibIdx{curLen+nn} = cat(2,neibIdx{curLen+nn}, bsIdxLen+tmpIdx);
                % neiFm{curLen+nn} = cat(1,neiFm{curLen+nn}, tmpIdx*0+j);
                Cij{curLen+nn} = cat(2,Cij{curLen+nn}, -log(1-D(nn,tmpIdx)/g.maxDist)*g.transitionFactor);% times distance factor
            end
        end
    end
    curLen = curLen+length(xCoord{i});
end


end