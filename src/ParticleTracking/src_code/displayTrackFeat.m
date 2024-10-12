function feature = displayTrackFeat(movieInfo, abSite, valLength, valTime, convergeTrackOnly)
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

ddCenter = abSite.loc.*abSite.resolution;%[286 346 12].*[0.415 0.415 1];
movieInfo.xCoord = movieInfo.xCoord*abSite.resolution(1);
movieInfo.yCoord = movieInfo.yCoord*abSite.resolution(2);
movieInfo.zCoord = movieInfo.zCoord*abSite.resolution(3);

%% redefine valid track
trackLength = cellfun(@length,movieInfo.tracks);
validTrack = find(trackLength>=valLength);

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
validTrack = validTrack(stopTimePt-stTime>valLength-1);
%% divide track into two classes
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
feature.valFlag = valFlag;
if convergeTrackOnly % remove those tracks that does not converge to ablation site after ablation
    validTrack(~valFlag) = [];
end
%% speed: first way consider jump as two separate adjacent edges
dist2ddCenter = zeros(length(validTrack),2);
confRatio = zeros(length(validTrack),3);
speed = zeros(length(validTrack),9); %before ablation speed, after, overall
allSpeed = nan(length(validTrack),stopTimePt-1);
%speedFrameWise = nan(length(validTrack), stopTimePt-1);
feature.val_frameIntervals = abSite.frameIntervals(1:stopTimePt-1);
for j=1:length(validTrack)
    trackNum = validTrack(j);
    curTrack = movieInfo.tracks{trackNum};
    curFrames = movieInfo.frames(curTrack);
    if curFrames(end) > stopTimePt
        curFrames = curFrames(curFrames<=stopTimePt);
        curTrack = curTrack(curFrames<=stopTimePt);
    end
    tmpDist2ddCenter = [movieInfo.xCoord(curTrack(1)) ...
        movieInfo.yCoord(curTrack(1)) movieInfo.zCoord(curTrack(1))];
    dist2ddCenter(j,1) = norm(tmpDist2ddCenter-ddCenter);
    tmpDist2ddCenter = [movieInfo.xCoord(curTrack(end)) ...
        movieInfo.yCoord(curTrack(end)) movieInfo.zCoord(curTrack(end))];
    dist2ddCenter(j,2) = norm(tmpDist2ddCenter-ddCenter);
    
    %confRatio(j,:) = confineRatio(movieInfo, trackNum);
    
    x = movieInfo.xCoord(curTrack);
    y = movieInfo.yCoord(curTrack);
    z = movieInfo.zCoord(curTrack);
    cc = [x y z];
    distPtWise = cc(2:end,:) - cc(1:end-1,:);
    distPtWise = sqrt(sum(distPtWise.^2,2));
    
    tmpT = zeros(length(curFrames)-1,1);
    for i=1:length(curFrames)-1
        tmpT(i) = sum(abSite.frameIntervals(curFrames(i):curFrames(i+1)-1));
    end
    spTmp = distPtWise./tmpT;
    for i=1:length(curFrames)-1
        allSpeed(j,curFrames(i):curFrames(i+1)-1) = spTmp(i);
    end
    
    if curFrames(1) < abSite.timePoint
        idx = find(curFrames <= abSite.timePoint);
        dist0 = sum(distPtWise(idx(1:end-1)));
        st = curFrames(idx(1));
        ed = curFrames(idx(end));
        t0 = sum(abSite.frameIntervals(st:ed-1));
    else
        dist0 = nan;
        t0 = nan;
    end
    if curFrames(end) > abSite.timePoint+1
        idx = find(curFrames >= abSite.timePoint+1);
        dist1 = sum(distPtWise(idx(1:end-1)));
        st = curFrames(idx(1));
        ed = curFrames(idx(end));
        t1 = sum(abSite.frameIntervals(st:ed-1));
    else
        dist1 = nan;
        t1 = nan;
    end

    distall = sum(distPtWise);
    tall = sum(tmpT);
    speed(j,:) = [dist0 dist1 distall t0 t1 tall dist0/t0 dist1/t1 distall / tall];%[dist0 / t0, dist1 / t1, distall / tall];
end
feature.dist2ddCenter = dist2ddCenter;
feature.speed = speed; % before, after, all
feature.confRatio = confRatio; % overall journay, overall distance, distance/journay
feature.abSite = abSite;
feature.allSpeed = allSpeed;
%% plot
% if 1 % do not display
%     dfEst = squeeze(nanmean(abs(vMatXYZ),1));
%     figure;plot(dfEst(:,1));
% 
%     vMatNorm = sqrt(nansum(vMatXYZ.^2,3));
%     vMatNorm(vMatNorm==0) = nan;
%     figure;plot(squeeze(nanmean(vMatNorm,1)))
%     xMat = vMatXYZ(:,:,1);
%     figure;plot(sum(~isnan(xMat),1));
% 
%     avgStraightV = sqrt(sum(avgStraightV.^2,2));
%     avgV = [];
%     for i=1:max(stTime)
%         avgV(i) = mean(avgStraightV(stTime==i));
%     end
%     figure;plot(avgV);
%     idx = stTime<30;
%     figure;scatter(dist2ddCenter(idx,1),confRatio(idx,2))
%     %% total valid particles
%     fms = cell(length(validTrack),1);
%     for j=1:length(validTrack)
%         trackNum = validTrack(j);
%         curTrack = movieInfo.tracks{trackNum};
%         fms{j} = movieInfo.frames(curTrack);
%     end
%     fmAll = cat(1, fms{:});
% end
end