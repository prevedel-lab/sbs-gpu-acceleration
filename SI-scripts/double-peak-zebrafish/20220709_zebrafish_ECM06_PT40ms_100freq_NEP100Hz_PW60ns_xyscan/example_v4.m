clear all;
% close all;
clc;

%% read the data 
%% Initialize variables.
data = dlmread('Data\PW60ns_x20y60z0um_xpt2ypt2z2um_PT40ms_100freq_NEP100Hz.txt','\t',5, 1);

% the piezo stage delay 40ms * 4 = 160 ms;
freq_num = 100;
AOM_diff = 0.04;
freq_begin = 4.65 - AOM_diff; % in GHz
freq_end = 7.5 - AOM_diff; % in GHz
pixel_X = 20*5+1;
pixel_Y = 60*5+1;
pixel_num = pixel_X * pixel_Y;
range_X = 20;
range_Y = 60;
X = linspace(0, range_X, pixel_X);
Y = linspace(0, range_Y, pixel_Y);

data_LIA = data(:,1);
data_LIA(freq_num*2000 : freq_num*2000 : length(data_LIA))= [];

% data_LIA(1 : freq_num*2000 : length(data_LIA))= [];

%% z plane number
z_plane_num = 1;
data_to_fit = data_LIA(pixel_X*pixel_Y*freq_num*(z_plane_num-1)+1:pixel_X*pixel_Y*freq_num*z_plane_num);

%%
fs = pixel_X*pixel_Y*freq_num;
% y = highpass(data_to_fit, freq_num*10, fs);
y = highpass(data_to_fit, 25, 25*100);
y_shift = 1;
data_to_fit = data_to_fit + y_shift;
y = y + y_shift;

figure;
plot(data_to_fit,'.')
hold on;
plot(y)
% 
% figure;
% pspectrum(data_to_fit,fs)
% hold on;
% pspectrum(y,fs)

%% To fit data
[data_X, data_Y, weights] = load_fixed_samples_per_pixel(y, freq_begin, freq_end, freq_num, pixel_num);

[parameters, states, chi_squares, number_iterations, execution_time] = lorentzian_gpufit_constrained(data_X, data_Y, freq_begin, freq_end, weights);

%%
figure;
% 'amplitude', 'shift', 'width' or 'offset'
displaying_fit("shift", data_X, data_Y, X, Y, weights, parameters, chi_squares, pixel_X, pixel_Y, true, z_plane_num)
