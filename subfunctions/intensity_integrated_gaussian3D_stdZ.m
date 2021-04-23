function intensity = intensity_integrated_gaussian3D_stdZ(xp,yp,zp,x,y,z,sxy,sz,dx,dy,dz)
    %computes intensity of pixels assuming a 3d gaussian PSF that is
    %integrated over each pixel; along z the envelope is that of a standard
    %gaussian (no integration).
    %(normalized Intensity, i.e. no prefactor)

    %xp,yp,zp, array of pixel indices
    %x,y,z: position of gaussian in real units (if pixels, enter dx=1 dy=1
    %dz=1)
    %convention: center (0,0,0) corresponds to corner of bottommost
    %z-plane
    %sxy: sigma xy in pixels
    %sz: sigma z in pixels
    %dx,dy,dz: voxel dimensions in real units
    
    %computing the center position of each voxel
    xp = ceil(xp) -0.5;
    yp = ceil(yp) -0.5;
    zp = round(zp);
    
    diffx1 =  (double(xp-0.5))*dx - x;
    diffx1 = diffx1 / ( sqrt(2.0)*sxy);% variable change for the erf function
    diffx2 =  (double(xp+0.5))*dx - x;
    diffx2 = diffx2 / ( sqrt(2.0)*sxy);

    diffy1 =  (double(yp-0.5))*dy - y;
    diffy1 = diffy1 / ( sqrt(2.0)*sxy );
    diffy2 =  (double(yp+0.5))*dy - y;
    diffy2 = diffy2 / ( sqrt(2.0)*sxy );

    gz = exp( - ( double(zp) *dz - z).^2 /(2*sz^2));
    
    intensity = ...
        abs( erf( diffx1) - erf(diffx2) ) .*...
        abs( erf( diffy1) - erf(diffy2) ) .*...
        gz;
end
