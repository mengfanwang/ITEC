function outIm = imdisplayWithROI3D(orgIm3d, roi3d)
% generate a 3D colorful data with ROI overlay
[h,w,z] = size(orgIm3d);
orgIm3d = scale_image(orgIm3d, 0,1);
outIm = zeros(h,w,3,z);
colorIntensity = 0.8;
cmap = colormap('jet');
cmap = cmap(randperm(size(cmap,1)),:);
clNum = size(cmap,1);
nParticle = max(roi3d(:));
for i=1:z
    %r =orgIm3d(:,:,i);
    g =orgIm3d(:,:,i);
    %b =orgIm3d(:,:,i);
    tmpROI = roi3d(:,:,i);
    r = zeros(h,w);
    b = zeros(h,w);
    for j=1:nParticle
        cCnt = max(1, floor(clNum*j/nParticle));
        if cmap(cCnt,1)+cmap(cCnt,3)<0.1
            cmap(cCnt,1) = 1;
        end
        r(tmpROI==j) = colorIntensity*cmap(cCnt,1);
        b(tmpROI==j) = colorIntensity*cmap(cCnt,3);
    end
    outIm(:,:,1,i) = r;
    outIm(:,:,2,i) = g;
    outIm(:,:,3,i) = b;
end
end