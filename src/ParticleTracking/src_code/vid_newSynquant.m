function [zMapRes, synIdRes] = vid_newSynquant(vid, minSz, maxSz)
% particle detection using synQuant for 3D video
% vid is YXZT video
if nargin < 2
    minSz = 20;
    maxSz = 200;
end
maxVid = max(vid(:));
if maxVid>255
    vid = round(255 * vid / maxVid);
end
[h,w,z,t] = size(vid);
zMapRes = zeros(h,w,z,t);
synIdRes = zeros(h,w,z,t);
q = paraQ3D(1, 0.8); % number of channels (always 1); ratio used for noise estimation, small value==>small variance
p = paraP3D(0.05, 0,0,255,minSz, maxSz, 0.2, 4);
vox_x = 2e-7; % useless indeed
for i=1:t
    fprintf('Start handling # %d/%d frame! \n', i, t);
    imgIn = vid(:,:,:,i);
    % detection for one channel
    [H1,W1,D1] = size(imgIn);
    datx = zeros(D1,H1*W1,'uint8');
    for ii=1:D1
        tmp = imgIn(:,:,ii)';
        datx(ii,:) = tmp(:);
    end
    det_res = ppsd3D(datx, W1, H1, vox_x, p,q);
    zMap = det_res.ppsd_main.zMap;
    zMap(zMap < 5) = 0;
    for zz = 1:D1
        zMapRes(:,:,zz,i) = zMap(zz,:,:);
    end
    synIdRes(:,:,:,i) = bwlabeln(zMapRes(:,:,:,i) > 0);
end