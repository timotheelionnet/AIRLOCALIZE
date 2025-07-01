function [Idog,dog] = DOGfilter3(I,s,w,padSize,kernSize)
    % Performs a difference of gaussian filter on an image or z-stack I
    % s: center wavelength of the filter
    % w: width of the filter
    
    %optional arguments (enter empty array to set to default)
    %padsize: size of padding (in pixels) added to the image; default = 30;
    %kernsize: size of convolution kernel; default = 21;
    
    % default size of the kernel set to encompass most of the main peak of the
    % kernel
    if isempty(kernSize)
        kernSize = round(6*s*w + 1); % three times the large gaussian sigma in both directions
        if mod(kernSize,2) == 0
            kernSize = kernSize + 1;
        end
    end
    
    % default padding ~ half the size of the kernel
    if isempty(padSize)
        padSize = floor(kernSize/4); 
    end
    
    g1 = fspecial('Gaussian', kernSize, max(0.5,s/w) );
    g2 = fspecial('Gaussian', kernSize, s*w );
    dog = g1 - g2;
    
    if ndims(I) ==3
        [nx,ny,~] = size(I); 
        B = padarray(I,[padSize,padSize,0],'symmetric');
        Idog = convn(double(B), dog, 'same');
        Idog = Idog(padSize+1:padSize+nx,padSize+1:padSize+ny,:);
    else
        [nx,ny] = size(I);
        B = padarray(I,[padSize,padSize],'symmetric');
        Idog = conv2(double(B), dog, 'same');
        Idog = Idog(padSize+1:padSize+nx,padSize+1:padSize+ny);
    end
end