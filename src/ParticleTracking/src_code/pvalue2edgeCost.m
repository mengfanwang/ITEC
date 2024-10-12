function edgeCost = pvalue2edgeCost(pValue, g)

if strcmp(g.costCalMethod, 'zscore')
    % way 1 pvalue to zscore; problem: there can be negative z-score
    edgeCost = -norminv(pValue);% norminv(1-p_oc/2);
elseif strcmp(g.costCalMethod, 'fisher')
    % way 2 pvalue to 2 degree of freedom chi-square
    edgeCost = -2*log(pValue);
else
    % way 3 pvalue to 1 degree of freedom chi-square:trick
    edgeCost = chi2inv(1-pValue,1);%
end

edgeCost = edgeCost+1e-6*rand(length(pValue),1);
end