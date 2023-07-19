function y = deriv_L1_L2(x, A1, g1, x1, A2, g2, x2)
% Derivative of double peak Lorentzian
% g1 and g2 are already squarred
    
    y = 2 .* A1 .* g1 .* (x1 - x) ./ (( g1 + (x1 -x).^2 ).^2) ...
        + 2 .* A2 .* g2 .* (x2 - x) ./ (( g2 + (x2 -x).^2 ).^2) ;

end