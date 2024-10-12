function curRatioXYZ = combinedVar(totalTps,jumpCoeff, tmpPtsPre)
% x6-y5+(x10-x6+y5-y1)/8 ~ N(0,sigma)
curRatioXYZ = zeros(3,1);
if ~isempty(tmpPtsPre) && isempty(find(tmpPtsPre==0, 1))
    % x-direction
    tmpX = 1 + jumpCoeff(1)/totalTps;
    curRatioXYZ(1,1) = tmpX^2 + (jumpCoeff(1)/totalTps)^2;
    % y-direction
    tmpY = 1 + jumpCoeff(2)/totalTps;
    curRatioXYZ(2,1) = tmpY^2 + (jumpCoeff(2)/totalTps)^2;
    % z-direction
    tmpZ = 1 + jumpCoeff(3)/totalTps;
    curRatioXYZ(3,1) = tmpZ^2 + (jumpCoeff(3)/totalTps)^2;
else
    curRatioXYZ(:,1) = 1;
end

end