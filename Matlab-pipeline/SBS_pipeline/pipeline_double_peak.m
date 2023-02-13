function [volume_single_peak, volume_double_peak, single_peak, double_peak] =...
    pipeline_double_peak(data_X, data_Y, experiment_settings, ...
    previous_peak_fit, constraints)
% - Pipeline double peak - 
% v0.3 : 09/09/22 
%   - added SNR
%   - added measure of fit difference 

    for z=1:size(data_X, 3)
        %% Constraints for the double peak fitting
        double_constraints = zeros([2*7, size(data_X, 2)], 'single');
        % Constraints ranges for L_1
        double_constraints(1, :) = constraints(1, :) ; %/2; % Min amplitude
        double_constraints(2, :) = constraints(2, :) ; %/2; % Max amplitude
        double_constraints(3, :) = constraints(3, :) ;%constraints(3, :); % Min x_0
        double_constraints(4, :) = constraints(4, :) ;%constraints(4, :); % Max x_0
        double_constraints(5, :) = constraints(5, :).^2 ; %/2; % Min gamma
        double_constraints(6, :) = constraints(6, :).^2 ; %/2; % Max gamma

        % Contraints ranges for L_2
        %temp = sort(previous_peak_fit(z).parameters(3,:));
        %delta_x0 = temp(round(0.90 * size(previous_peak_fit(z).parameters, 2))) ;
        delta_x0 = 0.5;
        double_constraints(7, :) = constraints(1, :) ; %/2; % Min amplitude
        double_constraints(8, :) = constraints(2, :) ; %/2; % Max amplitude
        double_constraints(9, :) = -delta_x0 ; %constraints(3, :) ;% -delta_x0 ;% Min x_0 
        double_constraints(10, :) = delta_x0 ; %constraints(4, :) ;% delta_x0 ;% Max x_0
        double_constraints(11, :) = constraints(5, :).^2 ; %/2; % Min gamma
        double_constraints(12, :) = constraints(6, :).^2 ; %/2; % Max gamma

        %Constraints ranges for offset
        double_constraints(13, :) = constraints(7, :); % Min offset
        double_constraints(14, :) = constraints(8, :); % Max offset        

        %% Determining which sample points are going to be used for the fit
        weights = ones(size(data_X)) ;
        weights(1: 1 + experiment_settings.samples_ignored(1) - 1,:) = 0 ;
        weights(end - experiment_settings.samples_ignored(2) + 1:end,:) = 0 ;

        double_weights = weights ;  
%         for j=1:size(data_X, 2)
%             double_weights(1+(j-1)*100 : j*100) = 10 * previous_peak_fit(z).parameters(1, j) ...
%                 .* previous_peak_fit(z).parameters(3, j)^2 ./ ...
%                 (previous_peak_fit(z).parameters(3,j)^2 + (data_X(:,j) - previous_peak_fit(z).parameters(2, j)).^2 ) + previous_peak_fit(z).parameters(4, j);
%         end
        
        peak_width = 2 ;
        %double_weights( data_X < (previous_peak_fit(z).parameters(2,:) - peak_width * previous_peak_fit(z).parameters(3,:))) = 0 ;
        %double_weights( data_X > (previous_peak_fit(z).parameters(2,:) + peak_width * previous_peak_fit(z).parameters(3,:))) = 0 ;
        %double_weights( data_X < 3.6) = 0 ;
        %double_weights( data_X > 5.8) = 0 ;

        %% Fitting 
        [single_peak(z).parameters, single_peak(z).states, single_peak(z).chi_squares, single_peak(z).number_iterations, single_peak(z).execution_time] = ...
            gpufit_lorentzian_constrained(data_X(:,:,z), data_Y(:,:,z), experiment_settings.freq_begin, experiment_settings.freq_end, double_weights, previous_peak_fit.parameters, constraints);

        % Either use abosulote or use relative lorentzians
        [double_peak(z).parameters, double_peak(z).states, double_peak(z).chi_squares, double_peak(z).number_iterations, double_peak(z).execution_time] = ...
            ... gpufit_double_lorentzian_constrained(data_X(:,:,z), data_Y(:,:,z), experiment_settings.freq_begin, experiment_settings.freq_end, double_weights, previous_peak_fit, double_constraints);
            gpufit_double_lorentzian_relative_constrained(data_X(:,:,z), data_Y(:,:,z), experiment_settings.freq_begin, experiment_settings.freq_end, double_weights, previous_peak_fit, double_constraints);

         %% Compute SNR
        noise_level=[];
        for n=1:size(data_X(1,:),2)
            low_lim = single_peak(z).parameters(2,n) - peak_width *single_peak(z).parameters(3,n);
            high_lim = single_peak(z).parameters(2,n) + peak_width *single_peak(z).parameters(3,n);
            %index = [find(data_X(:,n) > low_lim, 1) find(data_X(:,n) > high_lim, 1)];
            std_dev = std([data_Y(data_X(:,n) < low_lim ,n)' data_Y(data_X(:,n) > high_lim,n)']);
            %snr = [ snr (single_peak(z).parameters(1,n) / std_dev)] ;
            noise_level = [noise_level std_dev] ;
        end

        %% Converting into a 2d array
            % Double peak
        volume_double_peak.L1.amplitude(:,:,z)  = convert_to_2d(double_peak(z).parameters(1,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.L1.shift(:,:,z)      = convert_to_2d(double_peak(z).parameters(2,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.L1.width(:,:,z)      = convert_to_2d(double_peak(z).parameters(3,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.L2.amplitude(:,:,z)  = convert_to_2d(double_peak(z).parameters(4,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.L2.shift(:,:,z)      = convert_to_2d(double_peak(z).parameters(5,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.L2.width(:,:,z)      = convert_to_2d(double_peak(z).parameters(6,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.offset(:,:,z)        = convert_to_2d(double_peak(z).parameters(7,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume_double_peak.error(:,:,z)         = convert_to_2d(double_peak(z).chi_squares, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_double_peak.states(:,:,z)        = convert_to_2d(double_peak(z).states, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume_double_peak.time(z) = double_peak(z).execution_time;
        volume_double_peak.weights = double_weights;

            % Single peak
        volume_single_peak.amplitude(:, :, z)   = convert_to_2d(single_peak(z).parameters(1,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_single_peak.shift(:, :, z)       = convert_to_2d(single_peak(z).parameters(2,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume_single_peak.width(:, :, z)       = convert_to_2d(single_peak(z).parameters(3,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_single_peak.offset(:, :, z)      = convert_to_2d(single_peak(z).parameters(4,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume_single_peak.error(:,:,z)         = convert_to_2d(single_peak(z).chi_squares, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume_single_peak.states(:,:,z)        = convert_to_2d(single_peak(z).states, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume_single_peak.noise_level(:,:,z)   = convert_to_2d(noise_level, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume_single_peak.time(z) = single_peak(z).execution_time;

        index =  1:experiment_settings.pixel_X*experiment_settings.pixel_Y  ;
        volume_single_peak.raw_data.X(:, :, z)  = data_X ;
        volume_single_peak.raw_data.Y(:, :, z)  = data_Y ;
        volume_single_peak.raw_data.weights(:,:,z) = weights ;
        volume_single_peak.raw_data.index(:, :, z) = convert_to_2d(index, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;

        %% Measure of fit difference
        discrepancy = sum((previous_peak_fit(z).parameters - single_peak(z).parameters).^2,1) ;
        volume_single_peak.discrepancy(:,:,z) = convert_to_2d(discrepancy, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 


    end
    %%  Flipping every other plane
    % May only work with uneven number of pixels ( X and Y)
    volume_double_peak.L1.amplitude(:,:,1:2:end)  = rot90(volume_double_peak.L1.amplitude(:,:,1:2:end), 2) ; 
    volume_double_peak.L1.shift(:,:,1:2:end)      = rot90(volume_double_peak.L1.shift(:,:,1:2:end), 2) ; 
    volume_double_peak.L1.width(:,:,1:2:end)      = rot90(volume_double_peak.L1.width(:,:,1:2:end), 2) ; 
    volume_double_peak.L2.amplitude(:,:,1:2:end)  = rot90(volume_double_peak.L2.amplitude(:,:,1:2:end), 2) ; 
    volume_double_peak.L2.shift(:,:,1:2:end)      = rot90(volume_double_peak.L2.shift(:,:,1:2:end), 2) ; 
    volume_double_peak.L2.width(:,:,1:2:end)      = rot90(volume_double_peak.L2.width(:,:,1:2:end), 2) ; 
    volume_double_peak.offset(:,:,1:2:end)        = rot90(volume_double_peak.offset(:,:,1:2:end), 2) ;
    volume_double_peak.error(:,:,1:2:end)         = rot90(volume_double_peak.error(:,:,1:2:end), 2) ; 
    volume_double_peak.states(:,:,1:2:end)        = rot90(volume_double_peak.states(:,:,1:2:end), 2) ;
    
    volume_single_peak.amplitude(:, :, 1:2:end) = rot90(volume_single_peak.amplitude(:, :, 1:2:end), 2);
    volume_single_peak.shift(:, :, 1:2:end)     = rot90(volume_single_peak.shift(:, :, 1:2:end), 2);
    volume_single_peak.width(:, :, 1:2:end)     = rot90(volume_single_peak.width(:, :, 1:2:end), 2);
    volume_single_peak.offset(:, :, 1:2:end)    = rot90(volume_single_peak.offset(:, :, 1:2:end), 2);
    volume_single_peak.error(:, :, 1:2:end)     = rot90(volume_single_peak.error(:, :, 1:2:end), 2);
    volume_single_peak.states(:, :, 1:2:end)    = rot90(volume_single_peak.states(:, :, 1:2:end), 2);
    volume_single_peak.noise_level(:, :, 1:2:end)   = rot90(volume_single_peak.noise_level(:, :, 1:2:end), 2);
    volume_single_peak.discrepancy(:, :, 1:2:end)   = rot90(volume_single_peak.discrepancy(:, :, 1:2:end), 2);

    volume_single_peak.raw_data.index(:, :, 1:2:end) = rot90(volume_single_peak.raw_data.index(:, :, 1:2:end), 2);

    volume_single_peak.snr(:, :, :) = volume_single_peak.amplitude ./ volume_single_peak.noise_level ;
end