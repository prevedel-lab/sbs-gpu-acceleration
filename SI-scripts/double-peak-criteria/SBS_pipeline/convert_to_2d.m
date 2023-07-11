function array_2d = convert_to_2d(array_1d, pixel_X, pixel_Y, shift) 
    %array_2d = zeros(pixel_X + max(shift) * 2 +1, pixel_Y, 'single');
    array_2d = zeros(pixel_X + max(shift) * 2, pixel_Y, 'single');
    array_2d(shift+1:pixel_X+shift,:) = reshape(array_1d, pixel_X, pixel_Y) ;
    
    % flip y scan and shift the lines - need a theoretical formula ?
    for s=1:length(shift)
        array_2d(:,s:length(shift):end) = circshift(array_2d(:,s:length(shift):end),-shift(s)) ;
    end
    array_2d(:,2:2:end) = flipud(array_2d(:,2:2:end)) ;
    array_2d = array_2d.' ;
end 