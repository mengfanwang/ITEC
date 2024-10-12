function [uCur, curRatioSquare, jumpCoeff, devXYZSingle, E2VVarV] = predictPosAndVarRatio(curPt, tmpPtsPre, neiPt, tmpPtsPost, totalTps, movieInfo, velocity, g, stdEstFlag)
% predict the location of next location in neibor particle's frame
% stdEstFlag == 1 means we are estimating std
% stdEstFlag == 0 means we are estimating next location
curCoord = [movieInfo.xCoord(curPt) movieInfo.yCoord(curPt) movieInfo.zCoord(curPt)];
neiCoord = [movieInfo.xCoord(neiPt) movieInfo.yCoord(neiPt) movieInfo.zCoord(neiPt)];
E2VVarV = 0;
if g.timeJump
    if stdEstFlag
        if isfield(movieInfo, 'edgeJumpCoeff')
            neiIdx = movieInfo.nei{curPt}==neiPt;
            jumpCoeff = movieInfo.edgeJumpCoeff{curPt}(neiIdx,:);
        else
            k = movieInfo.frames(neiPt)-movieInfo.frames(curPt);
            jumpCoeff = [k k k];
        end
    else
        [jumpCoeff, E2VVarV] = jumpCoefficient1Edge(curPt, tmpPtsPre, neiPt, tmpPtsPost, totalTps, movieInfo, velocity);
    end
    uCur = curCoord + velocity.*jumpCoeff;
else
    jumpCoeff = nan;
    uCur = curCoord + velocity;%*(movieInfo.frames(neiPt) - movieInfo.frames(curPt));
end

muXYZ = neiCoord - uCur;
if g.stdCombined
    curRatioSquare = ones(3,2);
    if g.timeJump
        curRatioSquare(:,1) = combinedVar(totalTps, jumpCoeff, tmpPtsPre);
        curRatioSquare(:,2) = combinedVar(totalTps, jumpCoeff, tmpPtsPost);
    else
        tmpJumpCoeff = ones(1,3);
        curRatioSquare(:,1) = combinedVar(totalTps, tmpJumpCoeff, tmpPtsPre);
        curRatioSquare(:,2) = combinedVar(totalTps, tmpJumpCoeff, tmpPtsPost);
    end
    stdRatioEst2 = sqrt(sum(curRatioSquare,2));
    devXYZSingle = muXYZ./stdRatioEst2'; % make them following the same var
else
    curRatioSquare = ones(3,2);
    devXYZSingle = muXYZ;
end

