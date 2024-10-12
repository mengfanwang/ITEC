function InstanceSegmentation(data_path, tmp_path, q)
% addpath(genpath("./signalGenerator"));
% addpath(genpath("./newPrincipalCurvatureFunctions_Wei"));
% addpath(genpath("./detector_Wei"));
% addpath(genpath("./imageIO_Wei"));
% addpath(genpath("./boundaryRefine"));
% addpath(genpath("./SIPCv4"));

%% generator folder
if ~exist(fullfile(tmp_path, 'InstanceSeg_res_tif'),'dir')
    mkdir(fullfile(tmp_path, 'InstanceSeg_res_tif'));
end
if ~exist(fullfile(tmp_path, 'InstanceSeg_res'),'dir')
    mkdir(fullfile(tmp_path, 'InstanceSeg_res'));
end
if ~exist(fullfile(tmp_path, 'priCvt'),'dir')
    mkdir(fullfile(tmp_path, 'priCvt'));
end
%%
tif_files = dir(fullfile(data_path, '/*.tif'));
for i=1:numel(tif_files)
    %%
    fprintf('Instance segmentation processing %d/%d\n', i, numel(tif_files));
%     fprintf('Start load data! \n');
    tic;
    %%
    [~, ImName, ~] = fileparts(tif_files(i).name);
    filePath=fullfile(data_path, [ImName '.tif']);
    dat=tifread(filePath);
    %% get foreground
    FGPath=fullfile(tmp_path, 'synQuant_res', [ImName '.mat']);
    load(FGPath,"id_mat");
    id_mat(:,:,1) = 0;
    idMap=double(id_mat);
    BG=idMap==0;
%     BG=true(size(dat));
%     toc;
    %% get curvature
    [PC_raw,FGmap,curvature_all]=PrcplCrvtr_scaleInvariant_3D_v5(dat,q.smFactorLst,q.zRatio,BG);
%     toc;
    %% dilate the gap
    zThres=0;
    eig_res_3d=PC_raw-zThres;
    gapMap=imdilate(eig_res_3d>0,strel('disk',1));
    eig_res_3d(gapMap)=max(eig_res_3d(gapMap),1e-3);
    %% boundary refine seed option #1
    seedMap=getSeedMap(eig_res_3d,FGmap,q.maxSize);
    %% boundary refine
    [refine_res,threshold_res]=regionWiseAnalysis4d_Wei10(seedMap,eig_res_3d,dat);
    threshold_res=uint8(threshold_res/q.scaleTerm*255);
%     toc;
    %% save result
%     fprintf('Start save result! \n');
    out_dat=gray2RGB_HD(max(min(dat/q.scaleTerm,1),0));
    out_newId=label2RGB_HD(refine_res);
    out=out_dat*0.8+out_newId*0.5;
    tifPath=fullfile(tmp_path, 'InstanceSeg_res_tif', ImName);
    tifwrite(out, tifPath);

    priCvtPath=fullfile(tmp_path, 'priCvt', ImName);
    save(priCvtPath,"eig_res_3d","-v7.3");

    refine_res=uint32(refine_res);
    cellPath=fullfile(tmp_path, 'InstanceSeg_res', ImName);
    save(cellPath,"refine_res","threshold_res","-v7.3");
    toc;  
end
reset(gpuDevice());