function [edgeCost, stdEstFrom2Cur, pValue] = distance2EdgeCost(neiPosition, uCur,  curStd, neiStd, curRatioSquare,g, timeGap)
if nargin<7
    timeGap = 1;
end
% use chi-square
df = 0;
%chiSquTestStat = 0;
if g.stdCombined
    stdEstFrom2Cur = sqrt(curRatioSquare(:,1)'.*curStd.^2 + curRatioSquare(:,2)'.*neiStd.^2);
else
    stdEstFrom2Cur = curStd;
end
if 0 % way 1 strcmp(g.varEstMethod,'independent')
    if ~isempty(uNei) && isempty(find(isnan(neiStd) | neiStd==0, 1))
        chiSquTestStat = sum(((curPosition-uNei)./neiStd).^2) ;
        df = df+3;
    end
    if ~isempty(uCur)
        chiSquTestStat = chiSquTestStat + sum(((neiPosition-uCur)./curStd).^2);
        df = df+3;
    end
else
    muXYZ = neiPosition-uCur;
    % problem: variance use which one
    if ~isempty(muXYZ)
        chiSquTestStat = sum((muXYZ./stdEstFrom2Cur).^2, 2) ;
        df = df+3;
        pValue = 1-chi2cdf(chiSquTestStat, df);
        if g.directionWise % fisher method
            chiSquTestStat = (muXYZ./stdEstFrom2Cur).^2;
            pvDirectional = 1-chi2cdf(chiSquTestStat, 1);
            fisherTestSt = -2*(sum(log(pvDirectional),2)+log(pValue));
            pValue = 1-chi2cdf(fisherTestSt, 2*4);
        end
    else
        error('No muXYZ exist!!!\n');
    end
end
if isfield(g,'jumpCost')
    if length(g.jumpCost)>1
        pValue = pValue*g.jumpCost(timeGap);
        pValue = sum(min([g.jumpCost; g.jumpCost*0+pValue]));
        edgeCost = pvalue2edgeCost(pValue, g);
    elseif length(g.jumpCost)==1
        edgeCost = pvalue2edgeCost(pValue, g);
        edgeCost = edgeCost+(timeGap-1)*g.jumpCost;
    else
        edgeCost = pvalue2edgeCost(pValue, g);
    end
else
    edgeCost = pvalue2edgeCost(pValue, g);
end
end