function [uNei, uCur, curRatioSquare, velocity, jumpCoeff, E2VVarV] = estimateMeanFrom2Tracks(movieInfo, curPt, neiPt, g)
% uNei, mu estimated the mean from neighboring track
% uCur, ....from current track
%curPtPosinTrack, curDif, curVec, curTps
%% extract track information
% current track can be put outside of interation
curInfo = movieInfo.particle2track(curPt,:);
if ~isnan(curInfo(1))
    curTrack = movieInfo.tracks{curInfo(1)};
    curTps = movieInfo.frames(curTrack);
    curPtPosinTrack = curInfo(2);
    curVec = [movieInfo.xCoord(curTrack), movieInfo.yCoord(curTrack), movieInfo.zCoord(curTrack)];
else
    curTrack = curPt;
    curTps = movieInfo.frames(curTrack);
    curPtPosinTrack = 1;
    curVec = [movieInfo.xCoord(curTrack), movieInfo.yCoord(curTrack), movieInfo.zCoord(curTrack)];
end
% neighboring track
neiInfo = movieInfo.particle2track(neiPt,:);
if ~isnan(neiInfo(1))
    neiTrack = movieInfo.tracks{neiInfo(1)};
    neiTps = movieInfo.frames(neiTrack);
    neiPtPosinTrack = neiInfo(2);
    neixVec = movieInfo.xCoord(neiTrack);
    neiyVec = movieInfo.yCoord(neiTrack);
    neizVec = movieInfo.zCoord(neiTrack);
    neiVec = [neixVec neiyVec neizVec];
else
    neiTrack = neiPt;
    neiTps = movieInfo.frames(neiTrack);
    neiPtPosinTrack = 1;
    neiVec = [movieInfo.xCoord(neiTrack), movieInfo.yCoord(neiTrack), movieInfo.zCoord(neiTrack)];
end
[tmpPtsPre, tmpPtsPost] = determinePrePostPts(curPtPosinTrack, neiPtPosinTrack, neiTrack, g);


% way 2: estimate using average of all neighboring edges
uNei = [];
velocity = 0;
curRatioSquare = ones(3,2); % current variance ratio, neighbor variacne ratio in x,y,z directions
if isempty(tmpPtsPre) && isempty(tmpPtsPost)
    uCur = [];
    jumpCoeff = nan;
    E2VVarV = [];
else
    [velocity, ~, totalTps] = velocityCal(tmpPtsPre, tmpPtsPost, curVec,...
        curTps, curPtPosinTrack, neiVec, neiTps, neiPtPosinTrack);
    stdEstFlag = 0;% means we are estimating next location
    [uCur, curRatioSquare, jumpCoeff, ~, E2VVarV] = predictPosAndVarRatio(curPt, tmpPtsPre, ...
         neiPt, tmpPtsPost, totalTps, movieInfo, velocity, g, stdEstFlag);
     E2VVarV = [E2VVarV, length(tmpPtsPre)+length(tmpPtsPost)];
end

if isempty(uCur) && isempty(uNei)
    uCur = curVec(curPtPosinTrack,:);
end
end