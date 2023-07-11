function [parameters, states, chi_squares, number_iterations, execution_time] = gpufit_lorentzian_constrained(data_X, data_Y, freq_begin, freq_end, weights, initial_parameters, constraints)
%% Gpufit call
% Example of the Matlab binding of the Gpufit library implementing
% Levenberg Marquardt curve fitting in CUDA
% https://github.com/gpufit/Gpufit
%
% Simple example demonstrating a minimal call of all needed parameters for the Matlab interface
% http://gpufit.readthedocs.io/en/latest/bindings.html#matlab
tic
addpath("gpufit-matlab\")
if isempty(which('gpufit.m'))
    error('Gpufit library not found in Matlab path.');
end

assert(gpufit_cuda_available(), 'CUDA not available');

% data - convert to single precision
data_Y_single = single(data_Y);
data_X_single = single(data_X);
weights_single = single(weights);

% number of fits, number of points per fit
number_fits = size(data_X, 2);
number_samples = size(data_X, 1);

% model ID and number of parameter
model_id = ModelID.CAUCHY_LORENTZ_1D;
number_parameters = 4;

% initial parameters 
if (isempty(initial_parameters))
    initial_parameters = zeros(number_parameters, number_fits, 'single');
    for i=1:number_fits
        n_points = sum(weights(:,i)); 
        % Amplitude
        initial_parameters(1, i) = 1 ;

        % x_0
        initial_parameters(2, i) = (freq_begin + freq_end)/2-0.2;

        % gamma
        initial_parameters(3, i) = (freq_end - freq_begin) / 7.5;

        % offset
        initial_parameters(4, i) = 0.7;
    end
end

% Constrains

% amplitude
if (isempty(constraints))
    constraints = zeros([2*number_parameters, number_fits], 'single');
    constraints(1, :) = 0.35; % Min amplitude
    constraints(2, :) = 1.25; % Max amplitude
    constraints(3, :) = freq_begin + 0.2; % Min x_0
    constraints(4, :) = freq_end - 0.1; % Max x_0
    constraints(5, :) = 0.2; % Min gamma
    constraints(6, :) = 0.7; % Max gamma
    constraints(7, :) = 0.4; % Min offset
    constraints(8, :) = 1; % Max offset
end
constraint_types = int32([ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER]);

% run Gpufit and compute 1 fitted curve
max_number_iterations = 300 ;
precision = 1e-5 ;
gpufit_preparation = toc ;

[parameters, states, chi_squares, number_iterations, execution_time] = gpufit_constrained(data_Y_single, weights_single, ...
    model_id, initial_parameters, constraints, constraint_types, precision, max_number_iterations, [], EstimatorID.LSE, data_X_single);
end