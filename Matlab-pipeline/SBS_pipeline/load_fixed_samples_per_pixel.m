function [data_X, data_Y, weights] = ...
    load_fixed_samples_per_pixel(raw_data, freq_begin, freq_end, samples_per_pixel, n_pixels, samples_ignored)
    data_Y = reshape(raw_data, samples_per_pixel, n_pixels);
    data_X = repmat((linspace(freq_begin, freq_end, samples_per_pixel))',1 , n_pixels);
    weights = ones(samples_per_pixel, n_pixels);
    weights(1:1+samples_ignored(1)-1,:) = 0 ;
    weights(end - samples_ignored(2) + 1:end,:) = 0 ;
end