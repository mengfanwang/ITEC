function [neibIdx, Cij, edgeJumpCoeff] = transitCostInitial(xCoord,yCoord,zCoord,frames, g,q)
% we can learn initial variance and mean from data
%% using all the data to initialize transition cost
neibIdx = cell(g.particleNum, 1);
CX = cell(g.particleNum, 1);
CY = cell(g.particleNum, 1);
CZ = cell(g.particleNum, 1);
timepts = numel(xCoord);
curLen = 0;
for i=1:timepts
    %disp(i);
    for j=i+1:i+g.k
        if j>timepts
            break;
        end
        DX = abs(createDistanceMatrix(xCoord{i},xCoord{j}));
        DY = abs(createDistanceMatrix(yCoord{i},yCoord{j}));
        DZ = abs(createDistanceMatrix(zCoord{i},zCoord{j}));
        bsIdxLen = length(cat(1, xCoord{1:j-1}));
        for nn = 1:size(DX,1)
            tmpIdx = find(DX(nn,:) < g.maxDistXYZ(1) & DY(nn,:) < g.maxDistXYZ(2) & DZ(nn,:) < g.maxDistXYZ(3));
            if ~isempty(tmpIdx)
                neibIdx{curLen+nn} = cat(2,neibIdx{curLen+nn}, bsIdxLen+tmpIdx);
                CX{curLen+nn} = cat(2,CX{curLen+nn}, DX(nn,tmpIdx));
                CY{curLen+nn} = cat(2,CY{curLen+nn}, DY(nn,tmpIdx));
                CZ{curLen+nn} = cat(2,CZ{curLen+nn}, DZ(nn,tmpIdx));
            end
        end
    end
    curLen = curLen+length(xCoord{i});
end
validPts = ~cellfun(@isempty,CX);
nearestDX = cellfun(@min,CX(validPts));
nearestDY = cellfun(@min,CY(validPts));
nearestDZ = cellfun(@min,CZ(validPts));
if strcmp(g.varEstMethod,'median')
    stdFromAllPts = [nanmedian(abs(nearestDX-nanmean(nearestDX))), ...
        nanmedian(abs(nearestDY-nanmean(nearestDY))), ...
        nanmedian(abs(nearestDZ-nanmean(nearestDZ)))];
else
    stdFromAllPts = [std(nearestDX) std(nearestDY) std(nearestDZ)];
end
if g.stdCombined
    stdFromAllPts = stdFromAllPts/sqrt(2);
end
Cij = cell(g.particleNum, 1);
edgeJumpCoeff = cell(g.particleNum, 1);
for i=1:numel(CX)
    if isempty(CX{i})
        continue;
    end
%     % start calculate cost Cij
    neipos = [CX{i}',CY{i}',CZ{i}'];
    curpos = neipos*0;
    Cij{i} = distance2EdgeCost(neipos,curpos, stdFromAllPts, stdFromAllPts, ones(3,2),g);
    if g.timeJump
        edgeJumpCoeff{i} = zeros(size(neipos,1), 3);
        for j=1:size(neipos,1)
            edgeJumpCoeff{i}(j,:) = frames(neibIdx{i}(j))-frames(i);
        end
    else
        edgeJumpCoeff{i} = ones(size(neipos,1), 3);
    end
end
fprintf('finish trajectory intialization with purely distance!\n');

end