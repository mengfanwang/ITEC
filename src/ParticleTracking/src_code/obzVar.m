function [totalObzVar, vVar] = obzVar(curPt, tmpPtsPre, neiPt, tmpPtsPost, totalTps, movieInfo)
% estimate the overall observation variance (view all edges are the same)
bestK = movieInfo.frames(neiPt) - movieInfo.frames(curPt);

totalObzVar = [0 0 0];
vVar = [0 0 0];
curStd = movieInfo.particleStd(curPt,:);
neiStd = movieInfo.particleStd(neiPt,:);
if ~isempty(tmpPtsPre)
    totalObzVar = totalObzVar + curStd.^2*(1/totalTps+1/bestK)^2 + curStd.^2/totalTps^2;
    vVar = vVar+2*curStd.^2/totalTps^2;
else
    totalObzVar = totalObzVar + curStd.^2/totalTps^2;
end
if ~isempty(tmpPtsPost)
    totalObzVar = totalObzVar + neiStd.^2*(1/totalTps+1/bestK)^2 + neiStd.^2/totalTps^2;
    vVar = vVar+2*neiStd.^2/totalTps^2;
else
    totalObzVar = totalObzVar + neiStd.^2/totalTps^2;
end

% if sum(totalObzVar>10 & bestK==1 & totalTps==8)
%     aaa = 1;
% end
end