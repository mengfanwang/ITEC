function [kMat,E2VVarV] = jumpCoefficient1Edge(curPt, tmpPtsPre, neiPt, tmpPtsPost, totalTps, movieInfo, velocity)
if isempty(tmpPtsPre) && isempty(tmpPtsPost)
    kMat = nan;
    return;
end

bestK = movieInfo.frames(neiPt) - movieInfo.frames(curPt);

totalVar = [0 0 0];
if ~isempty(tmpPtsPre)
    curStd = movieInfo.particleStd(curPt,:);
    totalVar = totalVar+2*curStd.^2/totalTps^2;
end
if ~isempty(tmpPtsPost)
    neiStd = movieInfo.particleStd(neiPt,:);
    totalVar = totalVar+2*neiStd.^2/totalTps^2;
end
if isfield(movieInfo,'s')
    zz = movieInfo.s(:,1)==length(tmpPtsPre) & movieInfo.s(:,2)==length(tmpPtsPost);
    totalVar = totalVar + movieInfo.s(zz,[3,5,7]).*abs(velocity);
end
% we use velocity as E(velocity)
Ek2 = velocity.^2;

kMat = (bestK.*Ek2)./(Ek2+totalVar);
E2VVarV = [Ek2, totalVar];
end
% function [kMat, vObzVar] = jumpCoefficient1Edge(curPt, tmpPtsPre, neiPt, tmpPtsPost, totalTps, movieInfo, velocity)
% if isempty(tmpPtsPre) && isempty(tmpPtsPost)
%     kMat = nan;
%     return;
% end
% 
% bestK = movieInfo.frames(neiPt) - movieInfo.frames(curPt);
% totalObzVar = [0 0 0];
% vObzVar = [0 0 0];
% curStd = movieInfo.particleStd(curPt,:);
% neiStd = movieInfo.particleStd(neiPt,:);
% if ~isempty(tmpPtsPre)
%     totalObzVar = totalObzVar + curStd.^2*(1/totalTps+1/bestK)^2 + curStd.^2/totalTps^2;
%     vObzVar = vObzVar+2*curStd.^2/totalTps^2;
% else
%     totalObzVar = totalObzVar + curStd.^2/totalTps^2;
% end
% if ~isempty(tmpPtsPost)
%     totalObzVar = totalObzVar + neiStd.^2*(1/totalTps+1/bestK)^2 + neiStd.^2/totalTps^2;
%     vObzVar = vObzVar+2*neiStd.^2/totalTps^2;
% else
%     totalObzVar = totalObzVar + neiStd.^2/totalTps^2;
% end
% 
% preNum = length(tmpPtsPre);
% postNum = length(tmpPtsPost);
% tmpV2v1 = movieInfo.v2v1Var;
% idx = find(tmpV2v1(:,1)==preNum & tmpV2v1(:,2)==postNum & tmpV2v1(:,3)==bestK);
% if ~isempty(idx)
%     addVar = tmpV2v1(idx,4:6)-totalObzVar;
%     addVar(addVar<0)=0;
%     vObzVar = addVar+vObzVar;
% end
% % we use velocity as E(velocity)
% Ek2 = velocity.^2;
% 
% kMat = (bestK.*Ek2)./(Ek2+vObzVar);
% 
% end