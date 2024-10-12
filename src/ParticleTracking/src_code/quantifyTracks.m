function feature = quantifyTracks(movieInfo, abSite, valLength, valTime, convergeTrackOnly)
% extract the features from tracks

t = max(movieInfo.frames);
tmpTime = 0;
stopTimePt = t;
for i=abSite.timePoint+1:length(abSite.frameIntervals)
    tmpTime = tmpTime + abSite.frameIntervals(i);
    if valTime < tmpTime
        stopTimePt = i;
        break;
    end
end
feature.tpsBfAblation = abSite.timePoint;
feature.tpsAfAblation = stopTimePt - abSite.timePoint;
feature.ptsBfAblation = sum(movieInfo.frames<=abSite.timePoint) / feature.tpsBfAblation;
feature.ptsAfAblation = sum(movieInfo.frames<=stopTimePt & movieInfo.frames>abSite.timePoint) / feature.tpsAfAblation;

ddCenter = abSite.loc.*abSite.resolution;%[286 346 12].*[0.415 0.415 1];
movieInfo.xCoord = movieInfo.xCoord*abSite.resolution(1);
movieInfo.yCoord = movieInfo.yCoord*abSite.resolution(2);
movieInfo.zCoord = movieInfo.zCoord*abSite.resolution(3);

%% redefine valid track
trackLength = cellfun(@length,movieInfo.tracks);
validTrack = find(trackLength>=valLength);
if convergeTrackOnly % remove those tracks that does not converge to ablation site after ablation
    valFlag = ones(length(validTrack),1);
    for j=1:length(validTrack)
        trackNum = validTrack(j);
        curTrack = movieInfo.tracks{trackNum};
        curFrames = movieInfo.frames(curTrack);
        
        if curFrames(end) <= abSite.timePoint
            continue;
        else
            stPt = find(curFrames>=abSite.timePoint,1);
            endPt = find(curFrames>=stopTimePt,1);
            tmpDist2ddCenter = [movieInfo.xCoord(curTrack(stPt)) ...
                movieInfo.yCoord(curTrack(stPt)) movieInfo.zCoord(curTrack(stPt))];
            d1 = norm(tmpDist2ddCenter-ddCenter);
            tmpDist2ddCenter = [movieInfo.xCoord(curTrack(endPt)) ...
                movieInfo.yCoord(curTrack(endPt)) movieInfo.zCoord(curTrack(endPt))];
            d2 = norm(tmpDist2ddCenter-ddCenter);
            
            if d2 - d1 > 0
                valFlag(j) = 0;
            end
                
        end

    end
end
validTrack(~valFlag) = [];
%% # of tracks before and after ablation 
stTime = zeros(length(validTrack),1);
numTracks = zeros(length(abSite.frameIntervals)+1,1);
for j=1:length(validTrack)
    trackNum = validTrack(j);
    curTrack = movieInfo.tracks{trackNum};
    curFrames = movieInfo.frames(curTrack);
    
    stTime(j) = curFrames(1);
    numTracks(stTime(j):curFrames(end)) = numTracks(stTime(j):curFrames(end)) + 1;
end
validTrack = validTrack(stTime<=stopTimePt);
stTime = stTime(stTime<=stopTimePt);
feature.numTrackStartBfAblation = sum(stTime<=abSite.timePoint) / feature.tpsBfAblation;
feature.numTrackStartAfAblation = sum(stTime>abSite.timePoint) / feature.tpsAfAblation;
feature.numTrackBfAblation = mean(numTracks(1:abSite.timePoint));
feature.numTrackAfAblation = mean(numTracks(abSite.timePoint+1:stopTimePt));
feature.numTrackImmediateAfAblation = mean(numTracks(abSite.timePoint+1:...
    round(abSite.timePoint+feature.tpsAfAblation/5)));

%% speed: first way consider jump as two separate adjacent edges
%trackLength = trackLength(validTrack);
vMatXYZ = nan(length(validTrack), t-1, 3); % velocity for x,y,z direction
totalTime = zeros(length(validTrack),1);
avgStraightV = zeros(length(validTrack),3);
dist2ddCenter = zeros(length(validTrack),2);
confRatio = zeros(length(validTrack),3);
speed = zeros(length(validTrack),6); %before ablation speed, after, overall
for j=1:length(validTrack)
    trackNum = validTrack(j);
    curTrack = movieInfo.tracks{trackNum};
    curFrames = movieInfo.frames(curTrack);
    totalTime(j) = curFrames(end)-curFrames(1);
    tmpDist2ddCenter = [movieInfo.xCoord(curTrack(1)) ...
        movieInfo.yCoord(curTrack(1)) movieInfo.zCoord(curTrack(1))];
    dist2ddCenter(j,1) = norm(tmpDist2ddCenter-ddCenter);
    tmpDist2ddCenter = [movieInfo.xCoord(curTrack(end)) ...
        movieInfo.yCoord(curTrack(end)) movieInfo.zCoord(curTrack(end))];
    dist2ddCenter(j,2) = norm(tmpDist2ddCenter-ddCenter);
    
    avgStraightV(j,1) = (movieInfo.xCoord(curTrack(end))-...
        movieInfo.xCoord(curTrack(1)));
    avgStraightV(j,2) = (movieInfo.yCoord(curTrack(end))-...
        movieInfo.yCoord(curTrack(1)));
    avgStraightV(j,3) = (movieInfo.zCoord(curTrack(end))-...
        movieInfo.zCoord(curTrack(1)));
    
    confRatio(j,:) = confineRatio(movieInfo, trackNum);
    for k=1:length(curFrames)-1
        timeStart = curFrames(k);
        timeEnd = curFrames(k+1);
        gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
        ptStart = curTrack(k);
        ptEnd = curTrack(k+1);
        for tt = timeStart:timeEnd-1
            vMatXYZ(j, tt, 1) = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, tt, 2) = (movieInfo.yCoord(ptEnd)-movieInfo.yCoord(ptStart))/gapT;
            vMatXYZ(j, tt, 3) = (movieInfo.zCoord(ptEnd)-movieInfo.zCoord(ptStart))/gapT;
        end
    end
    timewiseDist = [movieInfo.xCoord(2:end)-movieInfo.xCoord(1:end-1)
        movieInfo.yCoord(2:end)-movieInfo.yCoord(1:end-1)
        movieInfo.zCoord(2:end)-movieInfo.zCoord(1:end-1)];
    timewiseDist = sqrt(sum(timewiseDist.^2, 2));
    
    dist0 = sum(timewiseDist(1:abSite.timePoint-1));
    t0 = sum(abSite.frameIntervals(1:abSite.timePoint-1));
    
    dist1 = sum(timewiseDist(abSite.timePoint+1:stopTimePt));
    t1 = sum(abSite.frameIntervals(abSite.timePoint+1:stopTimePt));
    
    distall = sum(timewiseDist(1:stopTimePt));
    tall = sum(abSite.frameIntervals(1:stopTimePt));
    speed(i,:) = [dist0, dist0 / t0, dist1, dist1 / t1, distall, distall / tall];
end
vMatNorm = sqrt(nansum(vMatXYZ.^2,3));
vMatNorm(vMatNorm==0) = nan;
avgSpeed = nanmean(vMatNorm,1);
feature.speedBfAblation = mean(avgSpeed(1:abSite.timePoint-1));
feature.speedAfAblation = mean(avgSpeed(abSite.timePoint+1:stopTimePt));
feature.speedImmediateAfAblation = mean(avgSpeed(abSite.timePoint+1:...
    round(abSite.timePoint+feature.tpsAfAblation/5)));
feature.speed = speed; % before, after, all
feature.confRatio = confRatio; % overall journay, overall distance, distance/journay
%% moving distances
minPtsNum = 0;
vMatXYZ = nan(length(validTrack), 2); % velocity for x,y,z direction
for j=1:length(validTrack)
    trackNum = validTrack(j);
    curTrack = movieInfo.tracks{trackNum};
    curFrames = movieInfo.frames(curTrack);
    tmpPos = find(curFrames>stopTimePt,1);
    if ~isempty(tmpPos)
        curTrack = curTrack(1:tmpPos-1);
        curFrames = curFrames(1:tmpPos-1);
    end
    if length(curTrack) < 2
        continue;
    end
    tmpPos = find(curFrames>abSite.timePoint,1);
    if isempty(tmpPos)
        timeStart = curFrames(1);
        timeEnd = curFrames(end);
        if timeEnd - timeStart > minPtsNum
            gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
            ptStart = curTrack(1);
            ptEnd = curTrack(end);

            tmpX = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpY = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpZ = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, 1) = sqrt(tmpX^2+tmpY^2+tmpZ^2);
        end
    elseif tmpPos == 1
        timeStart = curFrames(1);
        timeEnd = curFrames(end);
        if timeEnd - timeStart > minPtsNum
            gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
            ptStart = curTrack(1);
            ptEnd = curTrack(end);

            tmpX = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpY = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpZ = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, 2) = sqrt(tmpX^2+tmpY^2+tmpZ^2);
        end
    elseif tmpPos == length(curFrames)
        timeStart = curFrames(1);
        timeEnd = curFrames(end-1);
        if timeEnd - timeStart > minPtsNum
            gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
            ptStart = curTrack(1);
            ptEnd = curTrack(end-1);

            tmpX = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpY = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpZ = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, 1) = sqrt(tmpX^2+tmpY^2+tmpZ^2);
        end
    elseif tmpPos == 2
        timeStart = curFrames(2);
        timeEnd = curFrames(end);
        if timeEnd - timeStart > minPtsNum
            gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
            ptStart = curTrack(2);
            ptEnd = curTrack(end);

            tmpX = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpY = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpZ = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, 2) = sqrt(tmpX^2+tmpY^2+tmpZ^2);
        end
    else
        timeStart = curFrames(1);
        timeEnd = curFrames(tmpPos-1);
        if timeEnd - timeStart > minPtsNum
            gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
            ptStart = curTrack(1);
            ptEnd = curTrack(tmpPos-1);

            tmpX = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpY = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpZ = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, 1) = sqrt(tmpX^2+tmpY^2+tmpZ^2);
        end
        
        timeStart = curFrames(tmpPos);
        timeEnd = curFrames(end);
        if timeEnd - timeStart > minPtsNum
            gapT = sum(abSite.frameIntervals(timeStart:(timeEnd-1)));%timeEnd-timeStart;
            ptStart = curTrack(tmpPos);
            ptEnd = curTrack(end);

            tmpX = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpY = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            tmpZ = (movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart))/gapT;
            vMatXYZ(j, 2) = sqrt(tmpX^2+tmpY^2+tmpZ^2);
        end
    end
end
avgVelocity = nanmean(vMatXYZ,1);
feature.velocityBfAblation = avgVelocity(1);
feature.velocityAfAblation = avgVelocity(2);

%% second way only consider no jump adjacent
% t = max(movieInfo.frames);
% trackLength = cellfun(@length,movieInfo.tracks);
% validTrack = find(trackLength>=valLength); % should we?
% vMatXYZ = nan(length(validTrack), t-1, 3); % velocity for x,y,z direction
% 
% for j=1:length(validTrack)
%     trackNum = validTrack(j);
%     curTrack = movieInfo.tracks{trackNum};
%     curFrames = movieInfo.frames(curTrack);
%     adjLinks = find(curFrames(2:end)-curFrames(1:end-1)==1);
%     for k = 1:length(adjLinks)
%         timeStart = curFrames(adjLinks(k));
%         ptStart = curTrack(adjLinks(k));
%         ptEnd = curTrack(adjLinks(k)+1);
%         vMatXYZ(j, timeStart, 1) = movieInfo.xCoord(ptEnd)-movieInfo.xCoord(ptStart);
%         vMatXYZ(j, timeStart, 2) = movieInfo.yCoord(ptEnd)-movieInfo.yCoord(ptStart);
%         vMatXYZ(j, timeStart, 3) = movieInfo.zCoord(ptEnd)-movieInfo.zCoord(ptStart);
%     end
% end


%% plot
if 0 % do not display
    dfEst = squeeze(nanmean(abs(vMatXYZ),1));
    figure;plot(dfEst(:,1));

    vMatNorm = sqrt(nansum(vMatXYZ.^2,3));
    vMatNorm(vMatNorm==0) = nan;
    figure;plot(squeeze(nanmean(vMatNorm,1)))
    xMat = vMatXYZ(:,:,1);
    figure;plot(sum(~isnan(xMat),1));

    avgStraightV = sqrt(sum(avgStraightV.^2,2));
    avgV = [];
    for i=1:max(stTime)
        avgV(i) = mean(avgStraightV(stTime==i));
    end
    figure;plot(avgV);
    idx = stTime<30;
    figure;scatter(dist2ddCenter(idx,1),confRatio(idx,2))
    %% total valid particles
    fms = cell(length(validTrack),1);
    for j=1:length(validTrack)
        trackNum = validTrack(j);
        curTrack = movieInfo.tracks{trackNum};
        fms{j} = movieInfo.frames(curTrack);
    end
    fmAll = cat(1, fms{:});
end
end