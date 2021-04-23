function intensity = intensity_gaussian2D(xp,yp,x,y,sxy,dx,dy)
    %computes intensity of pixels assuming a standard 2d gaussian PSF 
    %(normalized Intensity, i.e. no prefactor)

    %xp,yp, array of pixel indices
    %x,y: position of gaussian center in real units (if pixels, enter dx=1 dy=1)
    %convention: center (0,0) corresponds to img corner 
    %sxy: sigma xy in pixels
    %dx,dy: pixel dimensions in real units
    
    %ceil(xp) - 0.5 is the coordinate of the center of the pixel the coordinate
    %xp belongs to
    xp = double(ceil(xp)) - 0.5;
    yp = double(ceil(yp)) - 0.5;
    
    gx = exp( - ( xp *dx - x).^2 /(2*sxy^2));
    gy = exp( - ( yp *dy - y).^2 /(2*sxy^2));
    intensity = gx.*gy;
end
