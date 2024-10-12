clc;clear;close all;
dbstop if error
%%% temporary code for visual something. Not in source

addpath('../dt');

addpath('../src_code_matlab');
addpath('../src_code_cellSegment');
addpath('../src_code_cellTracker');
addpath('../src_code_visualization'); 

%% system and path
load('/work/Mengfan/Embryo/TM0-49/res_reg_0.25/synQuant_refine_res_reg_wrong.mat');

for tt = 1:50
    tt_ind = num2str(99999+tt);
    tt_ind = tt_ind(2:6);
    data = tifread(['/work/Mengfan/Embryo/TM0-49/data_reg_0.25/' tt_ind '.tif']);
    threshold_reg{tt} = zeros(size(threshold_reg{tt}), 'uint8');
    for ii = 1:max(refine_reg{tt}(:))
        fprintf("%d:%d/%d\n", tt, ii, max(refine_reg{tt}(:)));
        label = find(refine_reg{tt} == ii);
        threshold_reg{tt}(label) = uint8(round(quantile(data(label),0.1)));
    end
end
save('/work/Mengfan/Embryo/TM0-49/res_reg_0.25/synQuant_refine_res_reg.mat', 'refine_reg','threshold_reg','-v7.3')