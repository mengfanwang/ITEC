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
movieInfo.orgCoord = [movieInfo.xCoord, movieInfo.yCoord, movieInfo.zCoord];
movieInfo.frames = cat(1, frames{:});
%% intial graph parameters
g = graphPara(length(movieInfo.xCoord));

movieInfo.Ci = zeros(g.particleNum,1)+g.observationCost; % use constant as observation cost
% initial transition cost
[neibIdx, Cij,edgeJumpCoeff] = transitCostInitial(xCoord,yCoord,zCoord,movieInfo.frames, g,q);
movieInfo.nei = neibIdx;
movieInfo.Cij = Cij;
movieInfo.edgeJumpCoeff = edgeJumpCoeff;

%% iterative update transition cost start
% the initial result from min-cost flow
movieInfoAll = cell(g.maxIter,1);

loopCnt = 1;
while 1
    disp(loopCnt);
    if loopCnt <=2 % result of first(or &second) iteration are initialization
        g.c_en = g.initEnter;% cost of appearance and disappearance in 1st/second run are 100
    else
        g.c_en = g.realEnter;% cost of appearance and disappearance in the scene
    end
    g.c_ex = g.c_en;
    g.observationCost = -(g.c_en+g.c_ex);
    
    % build graph using current transition cost
    [~, g, dat_in] = trackGraphBuilder(movieInfo, g);
    
    % min-cost circulation for min-cost flow
    movieInfo = mccTracker(dat_in, movieInfo, g);
    movieInfoAll{loopCnt} = movieInfo;
    
    % update jump Ratio==> punishment to jump
    g.jumpCost = movieInfo.jumpRatio;

    if loopCnt>g.maxIter
        break;
    end
    % there can be drift in the data, we should correct it from tracking results
    [movieInfo, dfEst] = driftFromTracks(movieInfo,g);
    if dfEst == 0 % no enough tracks
        break;
    end
    % update transition cost based on current trajectories
    [movieInfo, zz] = transitCostGaussinBidirection(movieInfo, g);

    loopCnt = loopCnt+1;
end

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
            %dataWithParticles = dataWithParticles(512:end, 512:end,:,:,:);
            save(fullfile(q.savePath,'dataWithParticles.mat'),'dataWithParticles','-v7.3');
        else
            load(fullfile(q.savePath,'dataWithParticles.mat'));
        end
    else
        if ~exist(fullfile(q.savePath,'dataWithParticlesOnly.mat'), 'file')
            dataWithParticles = generate3DImWithParticle(orgData, movieInfo, p);
            %dataWithParticles = dataWithParticles(512:end, 512:end,:,:,:);
            save(fullfile(q.savePath,'dataWithParticlesOnly.mat'),'dataWithParticles','-v7.3');
        else
            load(fullfile(q.savePath,'dataWithParticlesOnly.mat'));
        end
    end
    write3DVidwithTracks(orgData, movieInfo,  fullfile(q.savePath, q.trackFileName),q.orgDataName, p);
end
end