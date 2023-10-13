%% Clean up from previous file
clear all;
addpath("..\..\Matlab-pipeline\SBS_pipeline\")
addpath("..\..\Matlab-pipeline\gpufit-matlab\")
addpath("..\..\Matlab-pipeline\saveastiff_4.5\")
% close all;
clc;

%% read the data from the file
folder = './' ;
file_txt = 'PW60ns_x16ypt16z0um_xpt1ypt1z10um_PT80ms_100freq_NEP100Hz.txt' ;
file_zip = [file_txt(1:end-3) 'zip'] ; % replacing .txt by .zip
if isfile([folder file_txt])
    data = dlmread([folder file_txt],'\t',5, 1);
else 
    unzip([folder file_zip], folder) ;
    data = dlmread([folder file_txt],'\t',5, 1);
end

%%  Initialize variables - acquisition settings
% the piezo stage delay 40ms * 4 = 160 ms;
experiment_settings.samples_ignored = [10 0] ; % ignore first *5* samples, and last *10* samples
experiment_settings.shift           = [5 5] ;
experiment_settings.x_pixel_to_um   = 0.1 ;
experiment_settings.y_pixel_to_um   = 0.1 ;

experiment_settings.freq_num    = 100;
AOM_diff = 0.04+0.1;

experiment_settings.range_X     = 16; %in um
experiment_settings.range_Y     = 16; %in um
experiment_settings.n_planes    = 1 ;
experiment_settings.freq_begin  = 2.89 - AOM_diff; % in GHz
experiment_settings.freq_end    = 6.166 - AOM_diff; % in GHz

experiment_settings.pixel_X     = experiment_settings.range_X / experiment_settings.x_pixel_to_um + 1;
experiment_settings.pixel_Y     = experiment_settings.range_Y / experiment_settings.y_pixel_to_um + 1;
experiment_settings.pixel_num   = experiment_settings.pixel_X * experiment_settings.pixel_Y;

data_LIA = data(:,1);
data_LIA(experiment_settings.freq_num*2000 : experiment_settings.freq_num*2000 + 1: length(data_LIA))= [];

start = tic ;
%% Preprocessing
% Filter data
fs = experiment_settings.pixel_X*experiment_settings.pixel_Y*experiment_settings.freq_num;
% y = highpass(data_to_fit, freq_num*10, fs);
y = highpass(data_LIA, 25, 25*100);
y_shift = 1;
data_LIA = data_LIA + y_shift;
y = y + y_shift;

[data_X, data_Y] = pipeline_preprocess(y, experiment_settings.freq_begin, experiment_settings.freq_end, experiment_settings.freq_num, experiment_settings.pixel_num, experiment_settings.n_planes) ;
%% Manual correction for the first 2000 points 
% Some problem with the first data batch (a missing bit)
data_X(:, 1:2000, 1) = data_X(:, 1:2000, 1) - (experiment_settings.freq_end - experiment_settings.freq_begin) / (experiment_settings.freq_num - 1);

%% Single peak fitting
initial_parameters = zeros(4, size(data_X, 2), 'single'); % Fits are done plane by plane
initial_parameters(1, :) = 1 ; % Amplitude
initial_parameters(2, :) = (experiment_settings.freq_begin + experiment_settings.freq_end) /2-0.1;  % x_0
initial_parameters(3, :) = (experiment_settings.freq_end - experiment_settings.freq_begin) / 7.5; % gamma
initial_parameters(4, :) = 0.7; % offset


constraints = zeros([2*4, size(data_X, 2)], 'single'); % Fits are done plane by plane
constraints(1, :) = 0.3; % Min amplitude
constraints(2, :) = 2.0; % Max amplitude
constraints(3, :) = experiment_settings.freq_begin + 0.85; % Min x_0
constraints(4, :) = experiment_settings.freq_end - 0.6; % Max x_0
constraints(5, :) = 0.18; % Min gamma
constraints(6, :) = 0.75; % Max gamma
constraints(7, :) = 0.2; % Min offset
constraints(8, :) = 1; % Max offset

[volume, single_fit] = pipeline_single_peak(data_X, data_Y, ...
    experiment_settings, initial_parameters, constraints );

%% Double peak fitting
[volume_single_peak, volume_double_peak, single_peak, double_peak] =...
    pipeline_double_peak(data_X, data_Y, experiment_settings, ...
    single_fit, constraints);

%% Computing statistics
[AIC] = pipeline_analyze_fit(data_X, data_Y, experiment_settings, volume_single_peak, volume_double_peak);

%% Displaying 
experiment_settings.total_processing_time = toc(start) ;
% volume_display(experiment_settings, volume, [], []);
volume_display(experiment_settings, volume_single_peak, volume_double_peak, AIC);

%% Export the data - saving to .tiff file
clear options;
options.overwrite = true;
%saveastiff(volume_double_peak.L1.shift .* AIC_mask + volume.shift .* ( 1 - AIC_mask), 'tiff_export\L1_shift.tif', options);
saveastiff(volume.shift , 'tiff_export\XZ_shift.tif', options);
saveastiff(volume.width , 'tiff_export\XZ_width.tif', options);
saveastiff(volume.amplitude, 'tiff_export\XZ_amplitude.tif', options);
saveastiff(volume.offset, 'tiff_export\XZ_offset.tif', options);