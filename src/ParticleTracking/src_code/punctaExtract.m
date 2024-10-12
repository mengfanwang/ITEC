function [roiMask, orgData, xCoord,yCoord,zCoord, pvalue, frames] = punctaExtract(q)
% load the particle information from label image

if isfield(q,'synIdRes')
    roiMask = q.synIdRes;
else
    dat = load(fullfile(q.savePath, q.particleDataName),'synIdRes');
    roiMask = dat.synIdRes;
end
if isfield(q,'orgData')
    orgData = q.orgData;
else
    if isfield(q, 'valCh')
        vol=imreadBF(fullfile(q.fPath, q.orgDataName),[],[],q.valCh,1);
    else
        vol=imreadBF(fullfile(q.fPath, q.orgDataName),[],[],1,1);
    end
    orgData = squeeze(vol{1});
end
[h,w,z,t] = size(roiMask);

if q.cropFlag% crop a region out
    roiMask = roiMask(1:round(h/2),1:round(w/2):end,:,:);
    orgData = orgData(1:round(h/2),1:round(w/2):end,:,:);
end
%%%%%%%%%%%%%%%%%%%% 3D particle detection results illustration
if q.visualizePuncta 
    minIntensity = q.minIntensity;
    filePath = fullfile(q.savePath, 'punctaVisualize');
    fileName = q.punctaFileName;
    if isfield(q, 'maxT')
        write3DImWithROI(orgData(:,:,:,1:q.maxT), roiMask(:,:,:,1:q.maxT), minIntensity, filePath,fileName);
    else
        write3DImWithROI(orgData, roiMask, minIntensity, filePath,fileName);
    end
end
%%%%%%%%%%%%%%%%%
if ~exist(fullfile(q.savePath,['firstRun3D_',q.particleDataName,'.mat']),'file') || q.particleNew
    % for observation cost if we have
    if exist(fullfile(q.savePath, 'zScore.txt'),'file') && exist(fullfile(q.savePath, 'particleIndex.txt'),'file')
        zScore = importdata(fullfile(q.savePath, 'zScore.txt'),'\n');
        particleIdx = importdata(fullfile(q.savePath, 'particleIndex.txt'),'\n');
        %zScore = cell(numel(zScore),1);
        %particleIdx = cell(numel(particleIdx),1);
        for i=1:numel(zScore)
            zScore{i}  = str2num(zScore{i});
            particleIdx{i}  = str2num(particleIdx{i});
        end
    end
    xCoord = cell(t,1);
    yCoord = cell(t,1);
    zCoord = cell(t,1);
    pvalue = cell(t,1);
    %zScore = cell(t,1);
    frames = cell(t,1);
    for i=1:t
        %disp(i);
        featMapFinal = roiMask(:,:,:,i);
        if q.roiBW
            featMapFinal = bwlabeln(featMapFinal,6);
        end 
        orgImage = orgData(:,:,:,i);
        % position
        featPropFinal = regionprops(featMapFinal,orgImage,'Area','PixelIdxList','WeightedCentroid','MeanIntensity'); %'Extrema'
        particleSz = [featPropFinal(:).Area]';
        valiIdx = find(vertcat(featPropFinal(:,1).MeanIntensity)>q.minIntensity ...
            & particleSz>=q.minSz);
        if q.visualizePuncta 
            
        end
        nFeats = length(valiIdx);
        temp = vertcat(featPropFinal.WeightedCentroid);
        zCoord{i} = temp(valiIdx,3);
        yCoord{i} = temp(valiIdx,2);
        xCoord{i} = temp(valiIdx,1);
        % zscore--currently use the the same value
        if exist(fullfile(q.savePath, 'zScore.txt'),'file') 
            pvalue{i} = 1-normcdf(zScore{i}(valiIdx)',0,1);
        else
            pvalue{i} = 1e-16+zeros(nFeats,1);
        end
        frames{i} = i*ones(numel(valiIdx),1);
    end
    pvalue = cat(1, pvalue{:});
    save(fullfile(q.savePath,['firstRun3D_',q.particleDataName,'.mat']),'xCoord','yCoord', 'zCoord', 'pvalue','frames');
else
    load(fullfile(q.savePath,['firstRun3D_',q.particleDataName,'.mat']));
end


%% a temp function to remove overlapped particles
for cc=1:numel(xCoord)
    coord = [xCoord{cc}, yCoord{cc}, zCoord{cc}];
    for i=1:size(xCoord{cc},1)
        for j=i+1:size(xCoord{cc},1)
            dis = coord(i,1:3)-coord(j,1:3);
            tmp = sum(dis.^2);
            if tmp<1.5
                xCoord{cc}(j) = [];
                yCoord{cc}(j) = [];
                zCoord{cc}(j) = [];
            end
        end
    end
end
fprintf('finish extracting puncta from detection results!\n');


end