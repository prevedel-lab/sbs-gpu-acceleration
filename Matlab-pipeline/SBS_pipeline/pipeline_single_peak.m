function [volume, single_fit] = pipeline_single_peak(data_X, data_Y, ...
    experiment_settings, initial_parameters, constraints )

    % Determining with sample points are going to be used for the fit
    weights = ones(size(data_X)) ;
    weights(1:1+experiment_settings.samples_ignored(1)-1,:) = 0 ;
    weights(end - experiment_settings.samples_ignored(2) + 1:end,:) = 0 ;
       
    for z=1:size(data_X, 3)
        % Applying the lorentzian fit
        [single_fit(z).parameters, single_fit(z).states, single_fit(z).chi_squares, single_fit(z).number_iterations, single_fit(z).execution_time] = ...
            gpufit_lorentzian_constrained(data_X(:,:,z), data_Y(:,:,z), experiment_settings.freq_begin, experiment_settings.freq_end, weights, initial_parameters, constraints);

        volume.amplitude(:, :, z)   = convert_to_2d(single_fit(z).parameters(1,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume.shift(:, :, z)       = convert_to_2d(single_fit(z).parameters(2,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume.width(:, :, z)       = convert_to_2d(single_fit(z).parameters(3,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ; 
        volume.offset(:, :, z)      = convert_to_2d(single_fit(z).parameters(4,:), experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume.error(:,:,z)         = convert_to_2d(single_fit(z).chi_squares, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume.states(:,:,z)        = convert_to_2d(single_fit(z).states, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
        volume.time(z) = single_fit(z).execution_time;

        index =  1:experiment_settings.pixel_X*experiment_settings.pixel_Y  ;
        volume.raw_data.X(:, :, z)  = data_X ;
        volume.raw_data.Y(:, :, z)  = data_Y ;
        volume.raw_data.weights(:,:,z) = weights ;
        volume.raw_data.index(:, :, z) = convert_to_2d(index, experiment_settings.pixel_X, experiment_settings.pixel_Y, experiment_settings.shift) ;
    end
    
    % May only work with uneven number of pixels ( X and Y)
    volume.amplitude(:, :, 1:2:end) = rot90(volume.amplitude(:, :, 1:2:end), 2);
    volume.shift(:, :, 1:2:end)     = rot90(volume.shift(:, :, 1:2:end), 2);
    volume.width(:, :, 1:2:end)     = rot90(volume.width(:, :, 1:2:end), 2);
    volume.offset(:, :, 1:2:end)    = rot90(volume.offset(:, :, 1:2:end), 2);
    volume.error(:, :, 1:2:end)     = rot90(volume.error(:, :, 1:2:end), 2);
    volume.states(:, :, 1:2:end)    = rot90(volume.states(:, :, 1:2:end), 2);

    volume.raw_data.index(:, :, 1:2:end) = rot90(volume.raw_data.index(:, :, 1:2:end), 2);
end