clc;clear;close all;

%% set input output and temporary file location
data_path = '../data/tiff';                % files should in %5d format start from '00000'
result_path = '../result';
tmp_path = '../tmp';

%% add folders to path
addpath('./cc_ImHandle');
addpath('./src_code_synquant');
addpath('./src_code_cellSegment');
addpath('./src_code_visualization');

addpath('./src_code_matlab');
addpath('./src_code_cellTracker');
addpath(genpath('./CINDA'));
addpath(genpath('./ParticleTracking'));
addpath('./dt');
addpath('./TGMM_wrapper/');

%% set parameters

% general params
scale_term = 300;                        % rescale data intensity roughly between 0-255
background_intensity = 0;                % intensity threshold to clip background 
z_resolution = 5.86;                     % the ratio between x/y-resolution and z

% synQuant detection parameters
paras_synQuant.visualization = true;    % save deteciton images to ./tmp/synQuant_res_tif
paras_synQuant.clipping = background_intensity;
paras_synQuant.minIntensity = 25;       % intensity threshold for detection. Cannot be smaller than background_intensity
paras_synQuant.downSampleScale = 2;     % x/y downsample for acceleration. All following params work on downsampled data
paras_synQuant.sigma = [3 3 1];         % Gaussian smoothing factor
paras_synQuant.minSize = 200;           % minimium detection size      
paras_synQuant.maxSize = 30000;         % maximum detetion size
paras_synQuant.zThreshold = 2;          % hypothesis testing zscore threshold

% instance segmentation parameters
paras_instSeg.scaleTerm = scale_term;
paras_instSeg.zRatio = z_resolution;
paras_instSeg.maxSize = 10000;          % maximum size of an individual cell
paras_instSeg.smFactorLst=1.2.^(0:7)*2; % principle curvature smooth factor

% motion flow estimation parameters
paras_motFlow.sigma = 1;                % Gaussian pyramid smoothing sigma
paras_motFlow.downSampleScale = 2;      % x/y downsample for acceleration and GPU memory. Generally set to 2.
paras_motFlow.layerNum = 3;             % Gaussian pyramid layers. (x&y size/dsScale/2^layerNum) should be integers. Generally set to 5.
paras_motFlow.padSize = [15 15 5];      % should be longer than max motion distance
paras_motFlow.useSegRes = true;         % use segmentation results to help estimation. Not suggest to set false.

% tracking parameters
paras_tracking = initial_q(paras_motFlow.downSampleScale, true); % error correction settings. down sample scale different from flow estimation not tested.
paras_tracking.scaleTerm = scale_term;
paras_tracking.stLoc = [];             % start location [x y z] if crop data
paras_tracking.szCrop = [];            % crop size if crop data
tps_to_process = [];                   % specific time points to process. [] for all tps
graphPara_cell(nan);                   % graphical tracking settings. Not use here and for modify only.


%% main function -- preprocessing
SynQuantDetection(data_path, tmp_path, paras_synQuant);

VarianceMapEstimation(data_path, tmp_path, scale_term, background_intensity);

addpath(genpath('./src_code_priCvt'));
InstanceSegmentation(data_path, tmp_path, paras_instSeg);
rmpath(genpath('./src_code_priCvt'));

batch_size = MotionFlowEstimation_stable(data_path, tmp_path, paras_motFlow);

%% main function -- tracking
paras_tracking.grid_size = 2^paras_motFlow.layerNum;
paras_tracking.batch_size = [50 50 4]; %batch_size;
if isempty(tps_to_process)
    tif_files = dir(fullfile(data_path, '*.tif'));
    tps_to_process = 0:length(tif_files)-1;
end
Tracking(data_path, tmp_path, result_path, tps_to_process, paras_tracking);