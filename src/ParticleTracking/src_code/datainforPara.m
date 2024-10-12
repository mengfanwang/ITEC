function q = datainforPara(datapath, dataName, punctaDatapath, punctaDataName, savepath, ...
    appendixStr, minSize, minIntensity)

% load detection results and orginal data
q.fPath = datapath;
q.orgDataName = dataName; % name of images
q.particlePath = punctaDatapath;
q.particleDataName =  punctaDataName;
q.particleNew = false; %used for re-run the puncta extract procedure
q.roiBW = true; % if the detection results is binary or labeled index
q.minSz = minSize;
q.minIntensity = minIntensity;
q.savePthPostFix = appendixStr; % the appendix for saved resulting images
q.cropFlag = false; % for debug usage, true: crop small region of test data
q.visualizeTracks = true;  % true: visualize the tracks which are saved in 3D images 
q.visualizePuncta = true; % true: visualize detection results which are saved in 3D images
q.savePath = savepath;% save path
if ~exist(q.savePath,'file')
    mkdir(q.savePath);
end
if q.cropFlag% crop a region out
    q.savePath = [savepath,q.savePthPostFix];% save path: _In50_Out50
end

if q.visualizeTracks
    q.trackFileName = 'tracks';
end
if q.visualizePuncta
    q.punctaFileName = 'puncta';
end
end