function [data_X, data_Y] = pipeline_preprocess(raw_data, freq_begin, freq_end, samples_per_pixel, n_pixels, n_planes)
    data_Y = reshape(raw_data, samples_per_pixel, n_pixels, n_planes);
    data_X = repmat((linspace(freq_begin, freq_end, samples_per_pixel))',1 , n_pixels, n_planes);
end