function [roiMask, orgData, xCoord,yCoord,zCoord, pvalue, frames, q] = loadParticleInfo(q)
% load the particle information from label image
q.savePath = q.fPath;% save path
if q.cropFlag% crop a region out
    q.savePath = [q.fPath,q.savePthPostFix];% save path: _In50_Out50
end
if ~exist(q.savePath,'file')
    mkdir(q.savePath);
end

if 1%~exist([q.fPath,'\',q.particleDataName,'_orgDataInfo.mat'],'file') || q.particleNew
    ppsdMat = tifread([q.fPath, q.particleDataName]);
    [h,w,~] = size(ppsdMat);

    %zSlices = 21;
    roiMask = reshape(ppsdMat,h,w,q.zSlices,[]);
    orgData = tifread([q.fPath, q.orgDataName, '.tif']);
    orgData = reshape(orgData,h,w,q.zSlices,[]);    

    if q.cropFlag% crop a region out
%         roiMask = roiMask(95:95+255, 617:617+255,:,:);
%         orgData = orgData(95:95+255, 617:617+255,:,:);
        roiMask = roiMask(512:end,512:end,:,:);
        orgData = orgData(512:end,512:end,:,:);
    end
    %save([q.fPath,'\',q.particleDataName,'_orgDataInfo.mat'],'roiMask','orgData','-v7.3');
else % too slow
    load([q.fPath,'\',q.particleDataName,'_orgDataInfo.mat']);
end
%%%%%%%%%%%%%%%%%%%% 3D particle detection results illustration
minIntensity = 20;
filePath = [q.fPath, q.orgDataName];
fileName = q.orgDataName;
write3DImWithROI(orgData, roiMask, minIntensity, filePath,fileName)
%%%%%%%%%%%%%%%%%
tic;
t = size(roiMask,4);
q.timepts = t;
if ~exist([q.fPath,'\firstRun3D_',q.particleDataName,'.mat'],'file') || q.particleNew
    % for observation cost
    if exist([q.fPath, '\zScore.txt'],'file') && exist([q.fPath, '\particleIndex.txt'],'file')
        zScore = importdata([q.fPath, '\zScore.txt'],'\n');
        particleIdx = importdata([q.fPath, '\particleIndex.txt'],'\n');
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
        valiIdx = find(vertcat(featPropFinal(:,1).MeanIntensity)>30 & particleSz>=20);
        nFeats = length(valiIdx);
        temp = vertcat(featPropFinal.WeightedCentroid);
        zCoord{i} = temp(valiIdx,3);
        yCoord{i} = temp(valiIdx,2);
        xCoord{i} = temp(valiIdx,1);
        % zscore--currently use the the same value
        if exist([q.fPath, '\zScore.txt'],'file') 
            pvalue{i} = 1-normcdf(zScore{i}(valiIdx)',0,1);
        else
            pvalue{i} = 1e-16+zeros(nFeats,1);
        end
        frames{i} = i*ones(numel(valiIdx),1);
    end
    pvalue = cat(1, pvalue{:});
    save([q.fPath,'\firstRun3D_',q.particleDataName,'.mat'],'xCoord','yCoord', 'zCoord', 'pvalue','frames');
else
    load([q.fPath,'\firstRun3D_',q.particleDataName,'.mat']);
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
fprintf('finish loading!\n');
end