function outData = generate3DImWithParticle(orgData, movieInfo, p)
% write the 3D data into a given folder
[h,w,z,t] = size(orgData);
outData = zeros(h,w, 3,z,t);
% remove the ROIs with smaller intensity level
orgDataDim = p.orgDataDim;
if ~isfield(p, 'maxT')
    p.maxT = inf;
end
movieInfo.xCoord = movieInfo.orgCoord(:,1);
movieInfo.yCoord = movieInfo.orgCoord(:,2);
movieInfo.zCoord = movieInfo.orgCoord(:,3);
parfor i=1:min(t, p.maxT)
    disp(i);
    orgDataSingle = orgData(:,:,:,i)/orgDataDim;
    outIm = imdisplayWithParticle3D(orgDataSingle, movieInfo, i, p);
    outData(:,:,:,:, i) = outIm;
end
end