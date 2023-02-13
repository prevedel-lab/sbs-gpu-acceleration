function [fit_statistics] = pipeline_analyze_fit(data_X, data_Y, volume_single_peak, volume_double_peak)

    for z=1:size(data_Y, 3) % for every plane
        N = sum(volume_double_peak(z).weights > 0) ;
        for x=1:size(volume_double_peak(z).error, 1)
           for y=1:size(volume_double_peak(z).error, 2)
                index = volume_single_peak.raw_data.index(x, y);
                if(index > 0)
                    % AIC
                    fit_statistics.aic.double(x, y, z)     = calc_AIC(volume_double_peak(z).error(x, y), N(index), 7 +1);
                    fit_statistics.aic.single(x, y, z)     = calc_AIC(volume_single_peak(z).error(x, y), N(index), 4 +1);                

                    % BIC
                    fit_statistics.bic.double(x, y, z)     = calc_BIC(volume_double_peak(z).error(x, y), N(index), 7 +1);
                    fit_statistics.bic.single(x, y, z)     = calc_BIC(volume_single_peak(z).error(x, y), N(index), 4 +1);

                    % model 1 = restricted = single peak
                    % model 2 = unrestricted = double peak
                    % ((RSS1 - RSS2) / (p2 - p1) ) / (RSS2 / (n - p2))
                    % p1 < p2
                    fit_statistics.Ftest.val(x,y,z)        = ...
                        ((volume_single_peak(z).error(x, y) - volume_double_peak(z).error(x, y)) / (7 - 4)) ...
                        / (volume_double_peak(z).error(x, y) / ( N(index) - 7)) ;
                    fit_statistics.Ftest.P(x,y,z)  = fcdf(fit_statistics.Ftest.val(x,y,z), 7 - 4, ( N(index) - 7)); % p2-p1, n-p2
                end
           end
        end

    end
    
    % Computing AIC weights for each model
    fit_statistics.aic.delta.double = fit_statistics.aic.double - fit_statistics.aic.single ;
    fit_statistics.aic.delta.single = fit_statistics.aic.single - fit_statistics.aic.double ;
    fit_statistics.aic.weight.double = exp( -0.5 .* fit_statistics.aic.delta.double) ./ (1 + exp (-0.5 .* fit_statistics.aic.delta.double)) ;
    fit_statistics.aic.weight.single = exp( -0.5 .* fit_statistics.aic.delta.single) ./ (1 + exp (-0.5 .* fit_statistics.aic.delta.single)) ; 

    fit_statistics.aic.double_or_single = fit_statistics.aic.double < fit_statistics.aic.single ;
    fit_statistics.aic.double_likelyhood = fit_statistics.aic.weight.double ./ fit_statistics.aic.weight.single ;
    
    % Computing BIC delta and mask
    fit_statistics.bic.delta = fit_statistics.bic.single - fit_statistics.bic.double ;
    fit_statistics.bic.double_or_single = fit_statistics.bic.double < fit_statistics.bic.single ;
    


end