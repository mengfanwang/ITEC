% experimental results of sigma and mu of order statistics
% size of particle varies from 30-500
eleNum = 500;
expNum = 10000;
sigma = zeros(eleNum,eleNum);
mu = zeros(eleNum,eleNum);
totalData = randn(eleNum*2, expNum);
[stData, orMat] = sort(totalData,1, 'descend');
% [~, idMat] = sort(orMat,1, 'ascend');
% widthMap = [0:eleNum*2-1]*expNum;
% widthMap = repmat(widthMap, expNum,1);
% idMat = idMat+widthMap;
for i=10:500 % particle
    disp(i);
    tic;
    parfor j=10:500 % neighbor
        disp(j);
%         tmpMat = idMat(:, 1:i+j);
%         tmpIdMap = false(expNum, eleNum*2);
%         tmpIdMap(tmpMat(:)) = true;
%         tmp = stData(tmpIdMap);
        tmpOrMat = orMat<=(i+j);
        tmp = stData(tmpOrMat);
        tmp = reshape(tmp,i+j,[]);
        %tmp = sort(tmp,1, 'descend');
        muUp = mean(tmp(1:i,:));
        muDown = mean(tmp(i+1:i+j,:));
        delta = muUp-muDown;
        muNew = mean(delta);
        sigmaNew = std(delta);
        sigma(i,j) = sigmaNew;
        mu(i,j) = muNew;
    end
    toc;
end
save('sigmaMu.mat','sigma','mu');
csvwrite('sigma.txt',sigma);
csvwrite('mu.txt',mu);