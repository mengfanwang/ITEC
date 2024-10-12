function [zMapRes, synIdRes] = vid_synquant(vid, minSz, maxSz)
% particle detection using synQuant for 3D video
% vid is YXZT video
addpath(genpath('src_synquant\synquant_v2_paper\'));

if nargin < 2
    minSz = 20;
    maxSz = 400;
end
[h,w,z,t] = size(vid);
zMapRes = zeros(h,w,z,t);
synIdRes = zeros(h,w,z,t);
%binRes = false(h,w,z,t);
for i=1:t
    fprintf('Start handling # %d/%d frame! \n', i, t);
    imgIn = vid(:,:,:,i);
    xstd = noise.estNoiseTop(imgIn,'pair',0.2); %
    [zMap,res] = runSynQuant3DSingle(imgIn,xstd,minSz,maxSz);

    zMapRes(:,:,:,i) = zMap;
    synIdRes(:,:,:,i) = res;
    %binRes(:,:,:,i) = res>0;
end