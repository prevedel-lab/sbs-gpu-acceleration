function [parameters, states, chi_squares, number_iterations, execution_time] = ...
    gpufit_double_lorentzian_constrained(data_X, data_Y, freq_begin, freq_end, weights, single_fit, constraints)
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
model_id = ModelID.DOUBLE_LORENTZIAN_1D;
number_parameters = 7;

% initial parameters 
initial_parameters = zeros(number_parameters, number_fits, 'single');
for i=1:number_fits
    n_points = sum(weights(:,i)); 

    % Lorentzian #1
    % Amplitude_1
    initial_parameters(1, i) = single_fit.parameters(1,i);
    % x_0_1 : start at the left
    initial_parameters(2, i) = single_fit.parameters(2,i) - 0.5;
    % gamma_squarred_1
    initial_parameters(3, i) = single_fit.parameters(3,i)^2;

    % Lorentian #2 
    % Amplitude_2
    initial_parameters(1, i) = single_fit.parameters(1,i) ;
    % x_0_2 : start at the right
    initial_parameters(2, i) = single_fit.parameters(2,i) + 0.5;
    % gamma_squarred_2
    initial_parameters(3, i) = single_fit.parameters(3,i)^2;
   
    % offset
    initial_parameters(4, i) = 0.7;
end
% Constrains

if (isempty(constraints)) 
    % amplitude
    constraints = zeros([2*number_parameters, number_fits], 'single');
    % Constraints ranges for L_1
    constraints(1, :) = 0.1; % Min amplitude
    constraints(2, :) = 0.8; % Max amplitude
    constraints(3, :) = freq_begin + 0.6; % Min x_0
    constraints(4, :) = freq_end - 0.7; % Max x_0
    constraints(5, :) = 0.1^2; % Min gamma
    constraints(6, :) = 0.2^2; % Max gamma
    
    % Contraints ranges for L_2
    constraints(7, :)  = 0.1; % Min amplitude
    constraints(8, :)  = 0.5; % Max amplitude
    constraints(9, :)  = -0.2; % Min x_0
    constraints(10, :) = 0.2; % Max x_0
    constraints(11, :) = 0.05; % Min gamma
    constraints(12, :) = 0.15^2; % Max gamma
    
    %Constraints ranges for offset
    constraints(13, :) = 0.4; % Min offset
    constraints(14, :) = 1; % Max offset
end

%Setting contraints types
constraint_types = int32([ ...
    ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER, ... %L_1
    ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER, ConstraintType.LOWER_UPPER, ... %L_2
    ConstraintType.LOWER_UPPER]); %offset

% run Gpufit and compute 1 fitted curve
max_number_iterations = 500 ;
precision = 1e-5 ;
gpufit_preparation = toc ;

[parameters, states, chi_squares, number_iterations, execution_time] = gpufit_constrained(data_Y_single, weights_single, ...
    model_id, initial_parameters, constraints, constraint_types, precision, max_number_iterations, [], EstimatorID.LSE, data_X_single);

function [b, a] = swap(a, b)
end

for i=1:size(number_iterations, 2)
    if(parameters(2,i) > parameters(5,i))
        [parameters([1,2,3],i), parameters([4,5,6],i)] = swap(parameters([1,2,3],i), parameters([4,5,6],i));
    
    end
end
end