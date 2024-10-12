function [speed, gap2ablation] = speedExtract(movieInfo, abSite, valLength, valTime, convergeTrackOnly)
% get the speed defined by valTime
% INPUT:
% OUTPUT:
% speed: a matrix contains the information related to each track, where
% each row correspond to one track and for each column:
% 1. the distance traveled before ablation
% 2. the speed before ablation
% 3. the distance traveled ~2min after ablation (defined in valTime)
% 4. the speed ~2min after ablation
% 5. the distance traveled ~10min after abaltion
% 6...
% PS: 2min and 10min may vary depending on the values in valTime
% gap2ablation: a matrix contains the relationship between each track and
% the ablation site, which are:
% 1. if it starts before ablation or after ablation
% 2. the time gap between its starting and ablation
% 3. distances (x-, y-, z-, overall) to the ablation site at the time of ablation

% contact: ccwang@vt.edu, 01/20/2020

t = max(movieInfo.frames);
stopTimePt = t + valTime*0; 
for ii = 1:length(valTime)
    tmpTime = 0;
    for i=abSite.timePoint+1:length(abSite.frameIntervals)
        tmpTime = tmpTime + abSite.frameIntervals(i);
        if valTime(ii) < tmpTime
            stopTimePt(ii) = i;
            break;
        end
    end
end

movieInfo.xCoord = movieInfo.xCoord*abSite.resolution(1);
movieInfo.yCoord = movieInfo.yCoord*abSite.resolution(2);
movieInfo.zCoord = movieInfo.zCoord*abSite.resolution(3);
% refine tracks
validTrack = extractValTracks(movieInfo, valLength, abSite, stopTimePt, convergeTrackOnly);

%% speed: first way consider jump as two separate adjacent edges
speed = zeros(length(validTrack),length(valTime)*2+2); %before ablation speed, after, overall
for j=1:length(validTrack)
    trackNum = validTrack(j);
    curTrack = movieInfo.tracks{trackNum};
    curFrames = movieInfo.frames(curTrack);

    timewiseDist = [movieInfo.xCoord(curTrack(2:end))-movieInfo.xCoord(curTrack(1:end-1)),...
        movieInfo.yCoord(curTrack(2:end))-movieInfo.yCoord(curTrack(1:end-1)),...
        movieInfo.zCoord(curTrack(2:end))-movieInfo.zCoord(curTrack(1:end-1))];
    timewiseDist = sqrt(sum(timewiseDist.^2, 2));
    
    singleS = zeros(1, length(valTime)*2+2); % baseline, valTime1, valTime2.. overall
    for i=1:length(curFrames)-1
        if curFrames(i) < abSite.timePoint
            singleS(1) = singleS(1) + timewiseDist(i);
            singleS(2) = singleS(2) + ...
                sum(abSite.frameIntervals(curFrames(i):curFrames(i+1)-1));
        elseif curFrames(i) > abSite.timePoint
            for kk = 1:length(stopTimePt) % we have two stop time: 2min and 10 min
                if curFrames(i) < stopTimePt(kk)
                    singleS(1+kk*2) = singleS(1+kk*2) + timewiseDist(i);
                    singleS(2+kk*2) = singleS(2+kk*2) + ...
                        sum(abSite.frameIntervals(curFrames(i):curFrames(i+1)-1));
                end
            end
        end
    end
    singleS(2:2:end) = singleS(1:2:end)./singleS(2:2:end);
    speed(j,:) = singleS;
end

inval_vec = isnan(speed(:,2)) & isnan(speed(:,4)) & isnan(speed(:,6));
speed(inval_vec,:) = [];

%% gap2ablation
gap2ablation = zeros(length(validTrack),11); %before ablation speed, after, overall
ddCenter = abSite.loc.*abSite.resolution; % ablation center
for j=1:length(validTrack)
    trackNum = validTrack(j);
    curTrack = movieInfo.tracks{trackNum};
    curFrames = movieInfo.frames(curTrack);
    singleD = zeros(1,11);
    
    % temporally
    if curFrames(1) < abSite.timePoint
        singleD(1) = 1;
        singleD(2) = sum(abSite.frameIntervals(curFrames(1):abSite.timePoint-1));
    elseif curFrames(1) == abSite.timePoint
        singleD(1) = 0;
        singleD(2) = 0;
    else
        singleD(1) = -1;
        singleD(2) = sum(abSite.frameIntervals(abSite.timePoint:curFrames(1)-1));
    end
    % spatially
    st_loc = [movieInfo.xCoord(curTrack(1)), movieInfo.yCoord(curTrack(1)),...
        movieInfo.zCoord(curTrack(1))];
    end_loc = [movieInfo.xCoord(curTrack(end)), movieInfo.yCoord(curTrack(end)),...
        movieInfo.zCoord(curTrack(end))];
    
    singleD(3:5) = st_loc - ddCenter;
    singleD(6) = norm(st_loc - ddCenter);
    singleD(7:9) = end_loc - ddCenter;
    singleD(10) = norm(end_loc - ddCenter);
    
    ab_time = find(curFrames == abSite.timePoint, 1);
    if ~isempty(ab_time)
        ab_loc = [movieInfo.xCoord(curTrack(ab_time)), movieInfo.yCoord(curTrack(ab_time)),...
            movieInfo.zCoord(curTrack(ab_time))];
        singleD(11) = norm(ab_loc - ddCenter);
    end
    gap2ablation(j,:) = singleD;
end

gap2ablation(inval_vec,:) = [];
end