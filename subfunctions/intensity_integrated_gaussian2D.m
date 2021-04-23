function intensity = intensity_integrated_gaussian2D(xp,yp,x,y,sxy,dx,dy)
    %computes intensity of pixels assuming an integrated 2d gaussian PSF:
    %intensity is integrated over each pixel
    %(normalized Intensity, i.e. no prefactor)
    
    %xp,yp, array of pixel indices
    %x,y: position of gaussian in real units (if pixels, enter dx=1 dy=1)
    %convention: center (0,0) corresponds to img corner 
    %sxy: sigma xy in pixels
    %dx,dy: pixel dimensions in real units
    
    %computing the center position of each voxel
    xp = double(ceil(xp)) - 0.5;
    yp = double(ceil(yp)) - 0.5;
    
    diffx1 =  (double(xp-0.5))*dx - x;
    diffx1 = diffx1 / ( sqrt(2.0)*sxy);% variable change for the erf function
    diffx2 =  (double(xp+0.5))*dx - x;
    diffx2 = diffx2 / ( sqrt(2.0)*sxy);

    diffy1 =  (double(yp-0.5))*dy - y;
    diffy1 = diffy1 / ( sqrt(2.0)*sxy );
    diffy2 =  (double(yp+0.5))*dy - y;
    diffy2 = diffy2 / ( sqrt(2.0)*sxy );

    intensity = abs( erf( diffx1) - erf(diffx2) ).* abs( erf( diffy1) - erf(diffy2) ); 
    
end
