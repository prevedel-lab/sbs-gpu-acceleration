%% Clean up from previous file
clear all;
addpath("..\..\Matlab-pipeline\SBS_pipeline\")
addpath("..\..\Matlab-pipeline\gpufit-matlab\")
% close all;
clc;

% /!\ Not the same double peak fitting pipeline /!\
% Constraints changed to specific values


%% read the data from the file
folder = 'Data\220803_immersionoil_240umspacer_PT80ms_100freq_NEP100Hz_zscan03\' ;
file_txt = 'PW60ns_x0y0z25um_x0y0zpt01um_PT80ms_100freq_NEP100Hz_07.txt' ;
file_zip = [file_txt(1:end-3) 'zip'] ; % replacing .txt by .zip
if isfile([folder file_txt])
    data = dlmread([folder file_txt],'\t',5, 1);
else 
    unzip([folder file_zip], folder) ;
    data = dlmread([folder file_txt],'\t',5, 1);
end

%%  Initialize variables - acquisition settings
% the piezo stage delay 40ms * 4 = 160 ms;

experiment_settings.samples_ignored = [0 0] ; % ignore first *5* samples, and last *10* samples
experiment_settings.shift           = [0] ;
experiment_settings.x_pixel_to_um   = 0.1 ;
experiment_settings.y_pixel_to_um   = 0.01 ;
experiment_settings.range_X     = 0;
experiment_settings.range_Y     = 25;
experiment_settings.n_planes    = 0 + 1 ;

experiment_settings.freq_num    = 100;
AOM_diff = 0.04 + 0.1;
experiment_settings.freq_begin  = 5.23 - AOM_diff; % in GHz
experiment_settings.freq_end    = 9.15 - AOM_diff; % in GHz
experiment_settings.pixel_X     = experiment_settings.range_X / experiment_settings.x_pixel_to_um + 1;
experiment_settings.pixel_Y     = experiment_settings.range_Y / experiment_settings.y_pixel_to_um + 1;
experiment_settings.pixel_num   = experiment_settings.pixel_X * experiment_settings.pixel_Y;

data_LIA = data(:,1);
data_LIA(experiment_settings.freq_num*2000 : experiment_settings.freq_num*2000 : length(data_LIA))= [];

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

%% Removing a certain number of points at the start
%removing = 180; %180
%data_X(:, end-removing+1:end) = [] ;
%data_Y(:, end-removing+1:end) = [] ;
%experiment_settings.pixel_Y = experiment_settings.pixel_Y - removing ;

%% Single peak fitting
initial_parameters = zeros(4, size(data_X, 2), 'single'); % Fits are done plane by plane
initial_parameters(1, :) = 1 ; % Amplitude
initial_parameters(2, :) = 7.20;  % x_0
initial_parameters(3, :) = (experiment_settings.freq_end - experiment_settings.freq_begin) / 7.5; % gamma
initial_parameters(4, :) = 0.65; % offset


constraints = zeros([2*4, size(data_X, 2)], 'single'); % Fits are done plane by plane
constraints(1, :) = 0.15; % Min amplitude
constraints(2, :) = 1.75; % Max amplitude
constraints(3, :) = 7.1; % Min x_0
constraints(4, :) = 7.3; % Max x_0
constraints(5, :) = 0.16; % Min gamma
constraints(6, :) = 0.75; % Max gamma
constraints(7, :) = 0.4; % Min offset
constraints(8, :) = 1; % Max offset

[volume, single_fit] = pipeline_single_peak(data_X, data_Y, ...
    experiment_settings, initial_parameters, constraints );

%% Displaying 
experiment_settings.total_processing_time = toc(start) ;
%volume_display(experiment_settings, volume, [], []);


%% single peak resolution - Down
clf
hold on
windows_size = 40 ;
scale_up = 50 ;
slope = [770 1225];
ydata = volume(1).amplitude;

pixel_to_um= experiment_settings.y_pixel_to_um;
    
    range = slope(1):slope(2) ;
    %gaussEqn    = 'a*exp(-((x-b)/c)^2)+d' ;
    erfEqn      = 'a * erf( (x-b) / c ) + d' ;

    startPoints = [1.5 mean(slope) 10 0] ;
    f = fit(range', ydata(range), erfEqn, 'Start', startPoints);

    plot(f, range', ydata(range));
figure;
plot(ydata(range), 'k*', 'MarkerSize', 8)
hold on;
plot(f(770:1225), 'b', 'LineWidth', 3)


    res_in_um = f.c * 2 * sqrt(log(2)) * pixel_to_um 
    f

