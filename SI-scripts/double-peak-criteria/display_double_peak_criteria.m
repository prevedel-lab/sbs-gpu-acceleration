function display_double_peak_criteria(coord, volume_single_peak, volume_double_peak, spectral_resolution, n, delta )

i = coord(2) ;
j = coord(1) ;
raw_index = volume_single_peak.raw_data.index(i, j) ;
raw_X = volume_single_peak.raw_data.X(:,raw_index);
raw_Y = volume_single_peak.raw_data.Y(:,raw_index);
plot(raw_X, raw_Y, '.', 'DisplayName', 'Raw data') ;

single_Y = volume_single_peak.amplitude(i,j) .* volume_single_peak.width(i,j)^2 ./ (volume_single_peak.width(i,j)^2 + (volume_single_peak.shift(i,j) - raw_X) .^2) ...
    + volume_single_peak.offset(i,j);
double_Y = volume_double_peak.L1.amplitude(i,j) .* volume_double_peak.L1.width(i,j) ./ (volume_double_peak.L1.width(i,j) + (volume_double_peak.L1.shift(i,j) - raw_X) .^2) ...
    + volume_double_peak.L2.amplitude(i,j) .* volume_double_peak.L2.width(i,j) ./ (volume_double_peak.L2.width(i,j) + (volume_double_peak.L2.shift(i,j) - raw_X) .^2) ...
    + volume_double_peak.offset(i,j);
hold on
plot(raw_X, single_Y, 'DisplayName', 'Single peak fit') ;
plot(raw_X, double_Y, 'DisplayName', 'Double peak fit') ;


yyaxis right
x = volume_double_peak.L1.shift(i,j)-n*spectral_resolution:spectral_resolution:volume_double_peak.L2.shift(i,j)+n*spectral_resolution ;
deriv = deriv_L1_L2(x, ...
            volume_double_peak.L1.amplitude(i,j), volume_double_peak.L1.width(i,j), volume_double_peak.L1.shift(i,j),...
            volume_double_peak.L2.amplitude(i,j), volume_double_peak.L2.width(i,j), volume_double_peak.L2.shift(i,j)) ;


plot(x, (deriv),'--', 'DisplayName', 'Derivative') ;
plot(x, (deriv + delta),':', 'DisplayName', 'Derivative + delta') ;
plot(x, (deriv - delta),':', 'DisplayName', 'Derivative - delta') ;
yline(0)


xline(x(1), '--', 'DisplayName', 'L1 shift - n*spectral resolution')
xline(volume_double_peak.L1.shift(i,j), '-.', 'DisplayName', 'L1 shift')
xline(volume_double_peak.L2.shift(i,j), '-.', 'DisplayName', 'L2 shift')
xline(x(end), '--', 'DisplayName', 'L2 shift + n*spectral resolution')
legend

end