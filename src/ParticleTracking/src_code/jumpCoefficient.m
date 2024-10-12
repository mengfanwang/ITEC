function MeanMat = jumpCoefficient(movieInfo,g)
% if there is one jump, d = 2*v, but in our data, the coefficient k<2. We
% should use data to get the best estimation of k.
if ~isfield(movieInfo, 'jumpCoeff')
    movieInfo.jumpCoeff = zeros(g.k,3);
    for i=1:g.k
        movieInfo.jumpCoeff(i,:) = i;
    end
end
trackLength = cellfun(@length,movieInfo.tracks);
validIdx = find(trackLength>=g.trackLength4var);
gapT = cell(length(validIdx),1);
gapXYZ = cell(length(validIdx),1);
gapV = cell(length(validIdx),1);
for k=1:length(validIdx)
    trackNum = validIdx(k);
    curTrack = movieInfo.tracks{trackNum};
    fms = movieInfo.frames(curTrack);
    gatTtemp = cell(g.k,1);
    gapXYZtemp = cell(g.k,1);
    gapVtemp = cell(g.k,1);
    for kk = 1:g.k
        gatTtemp{kk} = fms(kk+1:end)-fms(1:end-kk);
        muXYZ = zeros(length(curTrack)-kk,3);
        vXYZ = zeros(length(curTrack)-kk,3);
        for j=1:length(curTrack)-kk
            if fms(j+kk)-fms(j)>g.k
                continue;
            end
            [~, ~, ~, velocity] = ...
                estimateMeanFrom2Tracks(movieInfo, curTrack(j), curTrack(j+kk), g);
            muXYZ(j,1) = movieInfo.xCoord(curTrack(j+kk))-movieInfo.xCoord(curTrack(j));
            muXYZ(j,2) = movieInfo.yCoord(curTrack(j+kk))-movieInfo.yCoord(curTrack(j));
            muXYZ(j,3) = movieInfo.zCoord(curTrack(j+kk))-movieInfo.zCoord(curTrack(j));
            vXYZ(j,:) = velocity;
        end
        gapXYZtemp{kk} = muXYZ;
        gapVtemp{kk} = vXYZ;

    end
    gapT{k} = cat(1,gatTtemp{:});
    gapXYZ{k} = cat(1,gapXYZtemp{:});
    gapV{k} = cat(1,gapVtemp{:});
end

gapT = cat(1,gapT{:});
gapXYZ = cat(1,gapXYZ{:});
gapV = cat(1,gapV{:});
MeanMat = zeros(g.k,3);
numSamples = zeros(3,3);
for i=1:g.k
    for j = 1:3
        valIdx = gapT==i;
        numSamples(i,j) = sum(valIdx);
        MeanMat(i,j) = gapV(valIdx,j)\gapXYZ(valIdx,j);
        if i<3
            k = i;
            figure;scatter(gapV(valIdx,j), gapXYZ(valIdx,j));hold on;
            x = min(gapV(valIdx,j)):0.1:max(gapV(valIdx,j));
            y = k*x;
            plot(x,y,'LineWidth',2);hold on;
            k=MeanMat(i,j);
            y = k*x;
            plot(x,y,'LineWidth',2);hold on;title(k);
        end
    end
end
%disp(numSamples);
end