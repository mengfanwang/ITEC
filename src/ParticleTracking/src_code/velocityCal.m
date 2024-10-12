function [velocity, totalDist, totalTps, absV] = velocityCal(tmpPtsPre, tmpPtsPost, curVec,curTps, curPtPosinTrack, neiVec, neiTps, neiPtPosinTrack)
% estimate the velocity by distance/time
totalDist = 0;
totalTps = 0;
absDist = 0;
if ~isempty(tmpPtsPre)
    if tmpPtsPre(end)+1 ~= curPtPosinTrack
        fprintf('the end of prePts is not current node!!\n');
    end
    totalDist = curVec(tmpPtsPre(end)+1,:)-curVec(tmpPtsPre(1),:);
    absDist = sum(curVec(tmpPtsPre+1,:)-curVec(tmpPtsPre,:), 1);
    totalTps = curTps(tmpPtsPre(end)+1)-curTps(tmpPtsPre(1));
end
if ~isempty(tmpPtsPost)
    if tmpPtsPost(1) ~= neiPtPosinTrack
        fprintf('the start of postPts is not current neighboring node!!\n');
    end
    absDist = absDist + sum(neiVec(tmpPtsPost+1,:)-neiVec(tmpPtsPost,:),1);
    totalDist = totalDist + neiVec(tmpPtsPost(end)+1,:)-neiVec(tmpPtsPost(1),:);
    totalTps = totalTps + neiTps(tmpPtsPost(end)+1)-neiTps(tmpPtsPost(1));
end
velocity = totalDist/totalTps;
absV = absDist/totalTps;
end