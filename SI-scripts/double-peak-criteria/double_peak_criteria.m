function criteria = double_peak_criteria(volume_double_peak, spectral_resolution, n, delta )

criteria = volume_double_peak.L1.shift ;
for i = 1:size(volume_double_peak.L1.shift, 1)
    for j= 1:size(volume_double_peak.L1.shift, 2)
        x = volume_double_peak.L1.shift(i,j)-n*spectral_resolution:spectral_resolution:volume_double_peak.L2.shift(i,j)+n*spectral_resolution ;

        deriv = deriv_L1_L2(x, ...
            volume_double_peak.L1.amplitude(i,j), volume_double_peak.L1.width(i,j), volume_double_peak.L1.shift(i,j),...
            volume_double_peak.L2.amplitude(i,j), volume_double_peak.L2.width(i,j), volume_double_peak.L2.shift(i,j)) ;
        
        list = [size(find(diff(sign(deriv))),2)];
        for t = -delta:0.05:delta
            list = [list size(find(diff(sign(deriv + t))),2)];
        end
        criteria(i,j) = max(list);
        
    end
end

end
