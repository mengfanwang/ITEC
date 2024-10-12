function p = initial_display_p(numTracks)

    % display paramters
    p.shiftScale = [10 10 5];
    p.orgDataDim = 1.5;
    p.maxInt = 220;%max(orgIm3d(:));
    p.particleSize = 1;
    p.particleCl = [1 1 0.1]*p.maxInt;
    p.otherParticleCl = [0 0 1]*p.maxInt;% color for all particles other than those in the track
    p.stoppedParticleCl = [1 0 1]*p.maxInt;
    p.lineWidth = 0.1;
    p.RegLineCl = [0 1 0]*p.maxInt;
    p.stoppeLineCl = [1 0 0]*p.maxInt;
    %p.fileName = q.orgDataName;
    p.cmap = hsv(numTracks);
    p.cmap = p.cmap(randperm(numTracks),:);
    %p.cmap = p.cmap(showOrder,:)*p.maxInt;
    % % % we replace all particles as blocks
    p.simulateFlag = 0;
    p.withOtherTrackFlag = 0;
    
    
end