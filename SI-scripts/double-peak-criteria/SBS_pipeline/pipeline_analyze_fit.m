function [AIC] = pipeline_analyze_fit(data_X, data_Y, volume_single_peak, volume_double_peak)
% v0.2.1 - 12/09
%   - Added AIC likelyhood 
%   - Added "+1" as parameter count (for variance of residual)
% v0.2 - 09/09/22
%   - Added "ratio" measurement

    for z=1:size(data_Y, 3) % for every plane
        N = sum(volume_double_peak(z).weights > 0) ;
        for x=1:size(volume_double_peak(z).error, 1)
           for y=1:size(volume_double_peak(z).error, 2)
                index = volume_single_peak.raw_data.index(x, y);
                if(index > 0)
                    % AIC
                    AIC.raw.double(x, y, z)     = calc_AIC(volume_double_peak(z).error(x, y), N(index), 7 +1);
                    AIC.raw.single(x, y, z)     = calc_AIC(volume_single_peak(z).error(x, y), N(index), 4 +1);                

                    % BIC
                    AIC.bic.double(x, y, z)     = calc_BIC(volume_double_peak(z).error(x, y), N(index), 7 +1);
                    AIC.bic.single(x, y, z)     = calc_BIC(volume_single_peak(z).error(x, y), N(index), 4 +1);

                    % model 1 = restricted = single peak
                    % model 2 = unrestricted = double peak
                    % ((RSS1 - RSS2) / (p2 - p1) ) / (RSS2 / (n - p2))
                    % p1 < p2
                    AIC.Ftest.val(x,y,z)        = ...
                        ((volume_single_peak(z).error(x, y) - volume_double_peak(z).error(x, y)) / (7 - 4)) ...
                        / (volume_double_peak(z).error(x, y) / ( N(index) - 7)) ;
                    AIC.Ftest.P(x,y,z)  = fcdf(AIC.Ftest.val(x,y,z), 7 - 4, ( N(index) - 7)); % p2-p1, n-p2
                end
           end
        end

    end
    
    AIC.delta_aic.double = AIC.raw.double - AIC.raw.single ;
    AIC.delta_aic.single = AIC.raw.single - AIC.raw.double ;
    AIC.weight.double = exp( -0.5 .* AIC.delta_aic.double) ./ (1 + exp (-0.5 .* AIC.delta_aic.double)) ;
    AIC.weight.single = exp( -0.5 .* AIC.delta_aic.single) ./ (1 + exp (-0.5 .* AIC.delta_aic.single)) ;    
    AIC.double_or_single = AIC.raw.double < AIC.raw.single ;
    AIC.double_likelyhood = AIC.weight.double ./ AIC.weight.single ;
    
    AIC.bic.delta = AIC.bic.single - AIC.bic.double ;
    AIC.bic.double_or_single = AIC.bic.double < AIC.bic.single ;
    
    AIC.mask = AIC.raw.double ;
    AIC.mask(AIC.weight.double < 0.5) = 0 ;
    AIC.mask(AIC.weight.double > 0.5) = 1 ;
    
    AIC.sanity = abs((volume_double_peak.L2.shift - volume_single_peak.shift) - (volume_single_peak.shift -  volume_double_peak.L1.shift)) ;

    AIC.ratio = min(volume_double_peak.L1.amplitude ./ volume_double_peak.L2.amplitude, volume_double_peak.L2.amplitude ./ volume_double_peak.L1.amplitude );

    AIC.shift_delta = abs(volume_double_peak.L1.shift(:,:,:) - volume_double_peak.L2.shift(:,:,:));
    %asymetry = abs(abs(volume_single_peak.shift(:,:,1) - volume_double_peak.L1.shift(:,:,1)) - abs(volume_single_peak.shift(:,:,1) - volume_double_peak.L2.shift(:,:,1)));
    index = (AIC.shift_delta > 0.17) ... % below 0.25-0.3 can't distinguish (simulation) ;
        .* (AIC.ratio > 0.2 ) ... %SNR ~20 => fitting noise = relative_strenght ~0.2 (?)
        .* AIC.double_or_single ;  % has to be a 'true' double peak
    AIC.filter = index ;

end