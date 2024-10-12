function [tmpPtsPre, tmpPtsPost] = determinePrePostPts(curPtPosinTrack, neiPtPosinTrack, neiTrack, g)
% given curPt and neiPt and their position in tracks
% determine the particles that should be used for velocity calculation
validPreNum = g.validPre;
validPostNum = g.validPost;
if curPtPosinTrack-validPreNum<1
    postNum = (validPreNum-curPtPosinTrack+1)+ validPostNum;
    tmpPtsPre = 1:curPtPosinTrack-1;
    tmpPtsPost = neiPtPosinTrack:min(neiPtPosinTrack+postNum-1, size(neiTrack,1)-1);
elseif neiPtPosinTrack+validPostNum-1>size(neiTrack,1)-1
    preNum = validPreNum + (neiPtPosinTrack+validPostNum-1 - (size(neiTrack,1)-1));
    tmpPtsPre = max(1,curPtPosinTrack-preNum):curPtPosinTrack-1;
    tmpPtsPost = neiPtPosinTrack:size(neiTrack,1)-1;
else
    tmpPtsPre = (curPtPosinTrack-validPreNum):curPtPosinTrack-1;
    tmpPtsPost = neiPtPosinTrack:(neiPtPosinTrack+validPostNum-1);
end
if isfield(g, 'maxEdgeNum')
    if length(tmpPtsPre)>g.maxEdgeNum
        tmpPtsPre = tmpPtsPre(length(tmpPtsPre)-g.maxEdgeNum+1:length(tmpPtsPre));
    end
    if length(tmpPtsPost)>g.maxEdgeNum
        tmpPtsPost = tmpPtsPost(1:g.maxEdgeNum);
    end
end
end