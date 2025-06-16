function smooth = smooth_image_and_subtract_background8(img,params,varargin)
    % applies a bandpass filter to the image img
    % if image is a stack, each slice is 2D-filtered independently.
    % params: airlocalizeParams object
        % version 7 of the function uses the properties
        % params.numDim: whether image is 2D or 3D
        % params.filterLo: low frequency cutoff in length in units of psf
            % Sigma
        % params.filterHi: hi frequency cutiff length in pixel units
        % params.psfSigma: lateral extension of the gaussian point spread function
     
    p = inputParser();
    
    
    if params.numDim == 3 || params.numDim == 2
        
        disp('  smoothing image ...');
        
        %factors 1.5 is a heuristic attempt to reproduce results
        %from the Fourier filter given the same parameters
        smooth = 1.5*DOGfilter(img,...
            params.filterHi, params.filterLo*params.psfSigma(1), [],[]);
    else
        disp('data type error: smoothing only possible on 2d images or 3d stacks');
        smooth = 0;
        return;
    end
    
    % subtract the background
    smooth = smooth - mean(smooth(:));
    
end

