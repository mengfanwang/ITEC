function [devXYZ, stdXYZ,validTrack, priorPara] = getStdFromTracks(movieInfo,g)
% get the standard deviation from existing tracks
xCoord = movieInfo.xCoord;
yCoord = movieInfo.yCoord;
zCoord = movieInfo.zCoord;
numTrajectories = numel(movieInfo.tracks);
validTrack = 0;
devXYZ = cell(numTrajectories,1);
stdXYZ = nan(numTrajectories,3);
if ~isfield(g,'dependencyFactor') %
    g.dependencyFactor = 1.4;
end
for i=1:numTrajectories
    curTrack = movieInfo.tracks{i};
    if length(curTrack)<g.trackLength4var % only use those trajectories >=g.trackLength4var
        continue;
    end
    validTrack = validTrack+1;
    curTps = movieInfo.frames(curTrack);
    xVec = xCoord(curTrack);
    yVec = yCoord(curTrack);
    zVec = zCoord(curTrack);
    Vec = [xVec, yVec, zVec];
    % way 1 of estimate variance
    devXYZ{i} = zeros(length(curTrack)-1, 3);
    for j=1:length(curTrack)-1
        [tmpPtsPre, tmpPtsPost] = determinePrePostPts(j, j+1, Vec, g);
        [velocity, ~, totalTps] = velocityCal(tmpPtsPre, tmpPtsPost, Vec,curTps, j, Vec, curTps, j+1);
        stdEstFlag = 1; %means we are estimating std
        [~, ~, ~, devXYZSingle] = predictPosAndVarRatio(curTrack(j), tmpPtsPre, ...
            curTrack(j+1), tmpPtsPost, totalTps, movieInfo, velocity, g, stdEstFlag);
        devXYZ{i}(j,:) = devXYZSingle;
    end
    % the way of estimate standard deviation
    if g.truncatedGaussian>0 % truncate the samples, because there can be outliers
        finalStd = truncateSamples(devXYZ{i},[], g);
        stdXYZ(i,:) = finalStd;
    else
        stdXYZ(i,:) = sqrt((sum(devXYZ{i}.^2./(size(devXYZ{i},1)))));% does not use -1
    end
end
% estimate a gamma distribution of precisions as a prior distribution and
% re-estimate all variance
% tau = 1/sigma^2
priorPara = nan(3,2); % gamma prior for 3 direction
if g.varPrior>0
    trackLength = cellfun(@length, movieInfo.tracks);
    [stLength, ord] = sort(trackLength, 'descend');
    chsenTracks = ord(1:min(g.varPrior, length(ord)));
    if strcmp(g.priorType, 'gamma')
        for i=1:3 % 3 directions
            tmpTau = 1./(stdXYZ(chsenTracks,i).^2); % precision
            phat = gamfit(tmpTau);
            priorPara(i,1) = phat(1); % alpha in direction i
            priorPara(i,2) = 1/phat(2); % beta in direction i
        end
    elseif strcmp(g.priorType, 'sInvX2')
        [df, tau2] = fitScaledInverX2(stdXYZ(chsenTracks,:).^2);
        priorPara = [df tau2];
        %elseif strcmp(g.priorType, 'weiAvg')
        %    priorPara = [mean(stdXYZ(chsenTracks,:).^2)', var(stdXYZ(chsenTracks,:).^2)'];
    elseif strcmp(g.priorType, 'weiAvg') || strcmp(g.priorType, 'Gauss')
        [muX,sigmaX] = normfit(stdXYZ(chsenTracks,1));
        [muY,sigmaY] = normfit(stdXYZ(chsenTracks,2));
        [muZ,sigmaZ] = normfit(stdXYZ(chsenTracks,3));
        %figure;scatter(trackLength,stdXYZ(:,1));title(num2str(sum(~isnan(stdXYZ(:,1)))))
        %figure;histogram(stdXYZ(chsenTracks,1),20);
        priorPara = [muX, sigmaX^2; muY, sigmaY^2; muZ, sigmaZ^2];
    else
        error('Undefined prior function type!\n');
    end
    for i=1:numTrajectories
        if trackLength(i)<g.trackLength4var % only use those trajectories >=g.trackLength4var
            continue;
        end
        % the way of estimate standard deviation
        if g.truncatedGaussian>0 % truncate the samples, because there can be outliers
            finalStd = truncateSamples(devXYZ{i}, priorPara, g);
            stdXYZ(i,:) = finalStd;
        else
            stdXYZ(i,:) = stdFromPrior(devXYZ{i}, priorPara, g);
        end
    end
end

end