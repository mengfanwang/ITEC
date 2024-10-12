function [zMap, synId, fMap] = m_Synquant4Embryo_Paramater(imgIn, q)
% cell detection using synQuant for 3D image
% INPUT:
% vid: is YXZ 3D image
% q: minIntensity
% minSz: the minimum size of valid region (default 100 for 3d, 20 for 2d)
% maxSz: the maximum size of valid region
% OUTPUT:
% zMap: the zscore map of detected regions
% synId: the id map of the detection regions
% fMap: foreground maps
% contact: ccwang@vt.edu mengfanw@vt.edu 02/04/2020     

minSz = q.minSize;                
maxSz = q.maxSize;      
z_threshold = q.zThreshold;
% the following two params are not important and don't need to modify
minfill = 0.1;
maxWHRatio = 10;  


fMap =  imgIn > q.minIntensity;
maxVid = max(imgIn(:));
if maxVid>255
    imgIn = round(255 * imgIn / maxVid);
end
q = paraQ3D(1,0,0.8); % number of channels (always 1); way to combine multiple channel (always 0);ratio used for noise estimation, small value==>small variance
% INPUTS of paraP3D
% 1. fdr threshold; 2. z-score threshold; 3. lowest intensity; 4. highest
% intensity; 5. minimum size; 6. largest size; 7. min fill; 8. max WHRatio
p = paraP3D(0, 0,0,255,minSz, maxSz, minfill, maxWHRatio);
vox_x = 2e-7; % useless indeed
% detection for one channel
[H1,W1,D1] = size(imgIn);
datx = zeros(D1,H1*W1,'uint8');
for ii=1:D1
    tmp = imgIn(:,:,ii)';
    datx(ii,:) = tmp(:);
end
det_res = ppsd3D(datx, W1, H1, vox_x, p,q);
zMap1 = det_res.ppsd_main.zMap;
synId1 = det_res.ppsd_main.kMap;

zMap = zeros(size(imgIn));
synId = zeros(size(imgIn));
for i=1:size(zMap,3)
    zMap(:,:,i) = zMap1(i,:,:);
    synId(:,:,i) = synId1(i,:,:);
end

synId(zMap<=z_threshold | ~fMap) = 0;
s = regionprops3(synId, {'VoxelIdxList'});
cnt = 0;
synId = zeros(size(imgIn));
for i=1:numel(s.VoxelIdxList)
    if length(s.VoxelIdxList{i}) >= minSz
        cnt = cnt + 1;
        synId(s.VoxelIdxList{i}) = cnt;
    end
end
end


