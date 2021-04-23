function intensity = intensity_gaussian3D(xp,yp,zp,x,y,z,sxy,sz,dx,dy,dz)
    %computes intensity of pixels assuming a standard 3d gaussian PSF 
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
    xp = ceil(xp) - 0.5;
    yp = ceil(yp) - 0.5;
    zp = round(zp);
    
    gx = exp( - ( double(xp) *dx - x).^2 /(2*sxy^2));
    gy = exp( - ( double(yp) *dy - y).^2 /(2*sxy^2));
    gz = exp( - ( double(zp) *dz - z).^2 /(2*sz^2));
    intensity = gx.*gy.*gz;
end