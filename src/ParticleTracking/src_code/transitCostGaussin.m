function movieInfo = transitCostGaussin(movieInfo, g)
% update transition cost using existing trajectories using precursors
% 
xCoord = movieInfo.xCoord;
yCoord = movieInfo.yCoord;
zCoord = movieInfo.zCoord;
numTrajectories = numel(movieInfo.tracks);
%Dij = movieInfo.Dij;% the distances of neighbors
validPreNum = g.validPre;
zz = zeros(1e7,1);
zzCnt = 1;
for i=1:numTrajectories
    curTrack = movieInfo.tracks{i};
    xVec = xCoord(curTrack);
    xDif = xVec(2:end)-xVec(1:end-1);
    yVec = yCoord(curTrack);
    yDif = yVec(2:end)-yVec(1:end-1);
    zVec = zCoord(curTrack);
    zDif = zVec(2:end)-zVec(1:end-1);
    xDif = [0; xDif];
    yDif = [0; yDif];
    zDif = [0; zDif];
    % way 1 of estimate variance
    devXYZ = zeros(length(curTrack)-1, 3);
    for j=1:length(curTrack)-1
        tmpPts = max(j-validPreNum+1, 1):j;
        % problem: if jump particle, should we times the length of time gap
        uX = mean(xDif(tmpPts))+xVec(j);
        uY = mean(yDif(tmpPts))+yVec(j);
        uZ = mean(zDif(tmpPts))+zVec(j);
        % update
        nxtNode = curTrack(j+1);
        devXYZ(j,:) = [xCoord(nxtNode)-uX,...
            yCoord(nxtNode)-uY,...
            zCoord(nxtNode)-uZ];
    end
    if size(devXYZ,1)>1
        stdXYZ = std(devXYZ);
    else
        continue;
    end
    % way 2 of estimate variance: assume point are static
    % stdXYZ = [std(xDif) std(yDif) std(zDif)];
    % stdXYZ = std(sqrt(xDif.^2+yDif.^2+zDif.^2));
    
    % start update cost Cij
    for j=1:length(curTrack)
        tmpPts = max(j-validPreNum+1, 1):j;
        % problem: if jump particle, should we times the length of time gap
        uX = mean(xDif(tmpPts))+xVec(j);
        uY = mean(yDif(tmpPts))+yVec(j);
        uZ = mean(zDif(tmpPts))+zVec(j);
        % update
        curNode = curTrack(j);
        neiUp = movieInfo.nei{curNode};
        for nn=1:length(neiUp)
            %tmpDistance = [xCoord(neiUp(nn))-xCoord(curNode),...
            %                yCoord(neiUp(nn))-yCoord(curNode),...
            %                    zCoord(neiUp(nn))-zCoord(curNode)];
            neiPosition = [xCoord(neiUp(nn)),...
                yCoord(neiUp(nn)),...
                zCoord(neiUp(nn))];
            %tmpDistance=norm(tmpDistance);
            movieInfo.Cij{curNode}(nn) = -log(mvnpdf(neiPosition,[uX,uY,uZ],stdXYZ));
            zz(zzCnt) = mvnpdf(neiPosition,[uX,uY,uZ],stdXYZ);
            zzCnt = zzCnt+1;
        end
    end
end
end