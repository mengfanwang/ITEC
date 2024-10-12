function dist = confineRatio(movieInfo, trackNum)
    curTrack = movieInfo.tracks{trackNum};
    x = movieInfo.xCoord(curTrack);
    y = movieInfo.yCoord(curTrack);
    z = movieInfo.zCoord(curTrack);
    cc = [x y z];
    distPtWise = cc(2:end,:) - cc(1:end-1,:);
    distPtWise = sum(sqrt(sum(distPtWise.^2,2)));
    distT = norm(cc(end,:) - cc(1,:));
    
    dist = [distPtWise, distT, distT/distPtWise];
end