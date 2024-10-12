function write3DVidwithTracks(orgData,  movieInfo,  filePath,fileName, p)
% write the 3D data with all tracks and particles into a given folder
% the particles on the same track are the same color
%
if ~exist(filePath, 'file')
    mkdir(filePath);
end
[~,~,~,t] = size(orgData);
if ~isfield(p, 'maxT')
    p.maxT = inf;
end
% remove the ROIs with smaller intensity level
for i=1:min(t, p.maxT)
    disp(i);
    orgDataSingle = orgData(:,:,:,i)/p.orgDataDim;
    % build the 3D data
    [h,w, z] = size(orgDataSingle);
    orgData3d = zeros(h,w,3,z);
    for j=1:z
        orgData3d(:,:,1,j) = orgDataSingle(:,:,j);
        orgData3d(:,:,2,j) = orgDataSingle(:,:,j);
        orgData3d(:,:,3,j) = orgDataSingle(:,:,j);
    end
    outIm = imdisplayWithTrack3D(orgData3d, movieInfo, i, p);
    tifwrite(outIm/255, fullfile(filePath, [fileName, '_', num2str(i)]));
end
end