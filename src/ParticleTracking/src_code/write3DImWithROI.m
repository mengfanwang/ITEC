function write3DImWithROI(orgData, roiData,minIntensity,  filePath,fileName)
% write the 3D data into a given folder
[~,~,~,t] = size(orgData);
if ~exist(filePath,'file')
    mkdir(filePath);
end

% remove the ROIs with smaller intensity level
parfor i=1:t
    disp(i);
    roiDataSingle = roiData(:,:,:,i);
    if max(roiDataSingle(:))<2
        roiDataSingle =  bwlabeln(roiDataSingle,6);
    end
    orgDataSingle = orgData(:,:,:,i);
    nParticle = max(roiDataSingle(:));
    for j=1:nParticle
         idx = find(roiDataSingle==j);
         [y,x,z] = ind2sub(size(roiDataSingle),idx);
         maxUYX = unique([y,x],'rows');
         %fillRatio = length(maxUYX)/prod(max(maxUYX)-min(maxUYX));
        avgIntensity = mean(orgDataSingle(roiDataSingle==j));
        if avgIntensity<minIntensity || length(y)<20 %|| fillRatio<0.5
            roiDataSingle(roiDataSingle==j) = 0;
        end
    end
    outIm = imdisplayWithROI3D(orgDataSingle, roiDataSingle);
    tifwrite(outIm, fullfile(filePath, [fileName, '_', num2str(i)]));
end
end