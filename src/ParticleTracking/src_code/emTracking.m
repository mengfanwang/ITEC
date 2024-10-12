function [movieInfo,movieInfoAll] = emTracking(q)
% we propose a iterative tracking framework

%% load detection results and orginal data
[~, orgData, xCoord,yCoord,zCoord, ~, frames] = punctaExtract(q);

%% the structure to save results
movieInfo=struct('orgCoord',[],'xCoord',[],'yCoord',[], 'zCoord', [], 'Ci', [], ...
    'frames', [], 'nei', [],'neiFra',[], 'Cij', [], 'tracks',[], 'jumpRatio',[]);

movieInfo.xCoord = cat(1, xCoord{:});
movieInfo.yCoord = cat(1, yCoord{:});
movieInfo.zCoord = cat(1, zCoord{:});
movieInfo.frames = cat(1, frames{:});

% main function
[movieInfo,movieInfoAll] = mcfTracking(movieInfo, xCoord,yCoord,zCoord);

if q.visualizeTracks
    % post-processing
    % recover all coordinate
    movieInfo.xCoord = movieInfo.orgCoord(:,1);
    movieInfo.yCoord = movieInfo.orgCoord(:,2);
    movieInfo.zCoord = movieInfo.orgCoord(:,3);
    
    trackLength = cellfun(@length, movieInfo.tracks);
    % TIME
    maxTime =  movieInfo.frames(cellfun(@max, movieInfo.tracks));
    minTime = movieInfo.frames(cellfun(@min, movieInfo.tracks));
    % DISTANCE
    xMax =  movieInfo.xCoord(cellfun(@max, movieInfo.tracks));
    xMin = movieInfo.xCoord(cellfun(@min, movieInfo.tracks));
    yMax =  movieInfo.yCoord(cellfun(@max, movieInfo.tracks));
    yMin = movieInfo.yCoord(cellfun(@min, movieInfo.tracks));
    zMax =  movieInfo.zCoord(cellfun(@max, movieInfo.tracks));
    zMin = movieInfo.zCoord(cellfun(@min, movieInfo.tracks));
    dist = sqrt(sum(([xMax*0.415, yMax*0.415, zMax]-[xMin*0.415, yMin*0.415, zMin]).^2,2));
    if isfield(q, 'valLen')
        movieInfo.tracks = movieInfo.tracks(trackLength>=q.valLen);
        movieInfo.pathCost = movieInfo.pathCost(trackLength>=q.valLen); % no use indeed
    else
        movieInfo.tracks = movieInfo.tracks(trackLength>=3);
        movieInfo.pathCost = movieInfo.pathCost(trackLength>=3); % no use indeed
    end

    
    numTracks = numel(movieInfo.tracks);
    rng(10); % set the same seed for reproduction
    showOrder = randperm(numTracks);
    %% display
    % display paramters
    p.shiftScale = [10 10 5];
    p.orgDataDim = 1.5;
    p.maxInt = 220;%max(orgIm3d(:));
    p.particleSize = 1;
    p.particleCl = [1 1 0.1]*p.maxInt;
    p.otherParticleCl = [0 0 1]*p.maxInt;% color for all particles other than those in the track
    p.stoppedParticleCl = [1 0 1]*p.maxInt;
    p.lineWidth = 0.3;
    p.RegLineCl = [0 1 0]*p.maxInt;
    p.stoppeLineCl = [1 0 0]*p.maxInt;
    p.fileName = q.orgDataName;
    p.cmap = hsv(numTracks);
    p.cmap = p.cmap(showOrder,:)*p.maxInt;
    % % % we replace all particles as blocks
    p.simulateFlag = 0;
    p.withOtherTrackFlag = 0;
    if isfield(q, 'maxT')
        p.maxT = q.maxT;
    else
        p.maxT = inf;
    end
    if ~p.simulateFlag
        if ~exist(fullfile(q.savePath,'dataWithParticles.mat'), 'file')
            dataWithParticles = generate3DImWithParticle(orgData, movieInfo, p);
%             %dataWithParticles = dataWithParticles(512:end, 512:end,:,:,:);
%             save(fullfile(q.savePath,'dataWithParticles.mat'),'dataWithParticles','-v7.3');
%         else
%             load(fullfile(q.savePath,'dataWithParticles.mat'));
        end
    else
        if ~exist(fullfile(q.savePath,'dataWithParticlesOnly.mat'), 'file')
            dataWithParticles = generate3DImWithParticle(orgData, movieInfo, p);
%             %dataWithParticles = dataWithParticles(512:end, 512:end,:,:,:);
%             save(fullfile(q.savePath,'dataWithParticlesOnly.mat'),'dataWithParticles','-v7.3');
%         else
%             load(fullfile(q.savePath,'dataWithParticlesOnly.mat'));
        end
    end
    write3DVidwithTracks(orgData, movieInfo,  fullfile(q.savePath, q.trackFileName),q.orgDataName, p);
end
end