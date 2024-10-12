function validTrack = extractValTracks(movieInfo, valLength, abSite, stopTimePt, convergeTrackOnly)
%% redefine valid track
ddCenter = abSite.loc.*abSite.resolution;%[286 346 12].*[0.415 0.415 1];

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
            endPt = find(curFrames>=stopTimePt(end),1);
            if isempty(endPt)
                endPt = length(curFrames);
            end
            if endPt == stPt
                continue;
            end
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
    validTrack(~valFlag) = [];

end

end