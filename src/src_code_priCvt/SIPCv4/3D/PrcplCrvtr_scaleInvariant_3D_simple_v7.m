function [curvature_all]=PrcplCrvtr_scaleInvariant_3D_simple_v7(dat,smFactorLst,zRatio)
%scale invariant principal curvature for 3D data, v3.1
% Detect 2D+3D curvature.
% This version is based on Noise normalization.
%
% input:    dat:                3D input data
%           smFactorLst:        vector of smooth factor
% output:   curvature_merged:   scale-invariant principal curvature
%
% 9/29/2022 by Wei Zheng

fprintf('Start curvature calculation! \n');
dat=single(dat);
%% initialize
N_smNum=length(smFactorLst);
curvature_all=zeros([size(dat) N_smNum],class(dat));
%%
[X,Y,Z]=size(dat);
Roverlap=ceil(max(smFactorLst)*3);
PatchSzX=floor(gpuDevice().AvailableMemory/4/6/(Y*Z))-2*Roverlap;
patchNum=ceil(X/PatchSzX);
%% iteratively get optimal curvature
for smCnt=1:N_smNum
    smFactor=[1 1 1/zRatio]*smFactorLst(smCnt);
    curvature_temp=zeros(size(dat),"single");
    for patchCnt=1:patchNum
        
        xStartRaw=(patchCnt-1)*PatchSzX+1;
        xEndRaw=min(patchCnt*PatchSzX,X);
        xStart=xStartRaw-Roverlap;
        xEnd=xEndRaw+Roverlap;
        if patchCnt==1
            xStart=xStartRaw;
        end
        if patchCnt==patchNum
            xEnd=xEndRaw;
        end

        temp=getCurvature_3D_v4d2(dat(xStart:xEnd,:,:),smFactor);
        curvature_temp(xStartRaw:xEndRaw,:,:)=temp(xStartRaw-xStart+1:end-(xEnd-xEndRaw),:,:);
    end
    
    curvature_all(:,:,:,smCnt)=curvature_temp;
end

end