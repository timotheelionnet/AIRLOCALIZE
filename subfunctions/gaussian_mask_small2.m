function [x0,y0,z0,N0,err0,dist,it] = gaussian_mask_small2(data,spotCenter,params)
% runs a 2D/3D gaussian mask algorithm to localize and quantify the intensity
% of a fluorescent spot
%data is the image/stack
%all units in pix with origin at the edge of pixel i.e. leftmost corner
%center corresponds to y = 0.5
%spot_ctr = [x y z] : the guess center coordinates (or [x,y] in 2D)
%params is an airlocalizeParams object holding the fit parameters
    %params.psfType : 'gaussian' or 'integratedGaussian' or
        %'integratedGaussianStdZ' (last option only in 3D)
    %params.psfSigma(1) : PSF width (in pixels)
    %params.psfSigma(2) : PSF height (in pixels - ignored in 2D)
    %params.fittedRegionSize : range (in PSF width units) of the region around the spot used for fitting. 
    %params.maxIterations is the maximum number of iterations of the equations allowed
        %before convergence
    %params.tol is the tolerance on the convergence (in lateral pixel dimension units)

%% parse input
xs = double(spotCenter(1)); 
ys = double(spotCenter(2)); 
sxy = double(params.psfSigma(1)); 
tol = double(params.tol);
maxcount = params.maxIterations;
cutSize = double(params.fittedRegionSize);
cutwidth = cutSize*sxy;
dx = 1;
dy = 1;
if ndims(data) == 3
    zs = double(spotCenter(3));
    sz = double(params.psfSigma(2));
    cutwidth(2) = cutSize*sz;
    dz = 1;
else
    z0 = 0;
end

%% initialize arays
x0 = zeros(1,maxcount,'double'); 
y0 = zeros(1,maxcount,'double'); 
N0 = zeros(1,maxcount,'double');

it = 1;
x0(1,it) = xs;
y0(1,it) = ys;
N0(1,it) = 0;

x = x0(1,it);
y = y0(1,it);
dist = zeros(1,maxcount,'double'); 

if ndims(data) == 3
    z0 = zeros(1,maxcount,'double');
    z0(1,it) = zs;
    z = z0(1,it);
end

%% compute boundaries of the ROI over which the mask is applied
if ndims(data) == 3
    [~,ROIlimits] = compute_ROI_boundaries2(...
        size(data),[x0(1,it),y0(1,it),z0(1,it)],cutwidth,0,'small');
    %voxel indices (integer)
    xp_min = ceil(ROIlimits(1,1)); 
    xp_max = ceil(ROIlimits(2,1));
    yp_min = ceil(ROIlimits(1,2)); 
    yp_max = ceil(ROIlimits(2,2));
    zp_min = round(ROIlimits(1,3)); 
    zp_max = round(ROIlimits(2,3));
    [xp,yp,zp] = ndgrid(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max);
    sdata = double(data(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max));
else
    [~,ROIlimits] = compute_ROI_boundaries2(...
        size(data),[x0(1,it),y0(1,it)],cutwidth,0,'small');
    %pixel indices (integer)
    xp_min = ceil(ROIlimits(1,1)); 
    xp_max = ceil(ROIlimits(2,1));
    yp_min = ceil(ROIlimits(1,2)); 
    yp_max = ceil(ROIlimits(2,2));
    [xp,yp] = ndgrid(xp_min:xp_max,yp_min:yp_max);
    sdata = double(data(xp_min:xp_max,yp_min:yp_max));
end

%% loop through iterations
it = 2;
tol = tol*(dx+dy)/2.0; 
tmp = tol+1;
while (it <= maxcount && tmp > tol)
    
    if ndims(data) == 3
        switch params.psfType
            case 'gaussian'
                intensity = intensity_gaussian3D(...
                    xp,yp,zp,x,y,z,sxy,sz,dx,dy,dz);               
            case 'integratedGaussian' 
                intensity = intensity_integrated_gaussian3D(...
                    xp,yp,zp,x,y,z,sxy,sz,dx,dy,dz);
            case 'integratedGaussianStdZ' 
                intensity = intensity_integrated_gaussian3D_stdZ(...
                    xp,yp,zp,x,y,z,sxy,sz,dx,dy,dz);                
        end
        intsum = intensity .* sdata;
        intsum = sum(intsum(:));
        
        sumsum = intensity.*intensity;
        sumsum = sum(sumsum(:));
        
        sumx = (double(xp-0.5)).*intensity.*sdata;
        sumx = sum(sumx(:));
        
        sumy = (double(yp-0.5)).*intensity.*sdata;
        sumy = sum(sumy(:));
        
        sumz = double(zp).*intensity.*sdata;
        sumz = sum(sumz(:));
        
        if intsum <= 0 || sumsum == 0
            [x0(1,it),y0(1,it),z0(1,it),N0(1,it)] = deal(-1,-1,-1,-1);
        else
            x0(1,it) = double(sumx) / double(intsum);
            y0(1,it) = double(sumy) / double(intsum);
            z0(1,it) = double(sumz) / double(intsum);
            N0(1,it) = double(intsum) / double(sumsum); 
            
            if   location_out_of_ROI(...
                    [x0(1,it),y0(1,it),z0(1,it)],...
                    [xp_min,xp_max,yp_min,yp_max,zp_min,zp_max],[dx,dy,dz])
                   [x0(1,it),y0(1,it),z0(1,it),N0(1,it)] = ...
                       deal(-1,-1,-1,-1);
            end      
        end

        dist(1,it) = ...
            sqrt( (x-x0(1,it))^2 + (y-y0(1,it))^2 + (z-z0(1,it))^2 );
        x = x0(1,it); 
        y = y0(1,it); 
        z = z0(1,it);
        tmp = dist(1,it);
        if x0(1,it) == -1
            tmp = tol-1;
        end
        it = it+1;
        
    elseif ndims(data) == 2
        switch params.psfType
            case 'gaussian'
                intensity = intensity_gaussian2D(...
                    xp,yp,x,y,sxy,dx,dy);
            case 'integratedGaussian' 
                intensity = intensity_integrated_gaussian2D(...
                    xp,yp,x,y,sxy,dx,dy);
            case 'integratedGaussianStdZ' % option only available in 3D
                intensity = intensity_integrated_gaussian2D(...
                    xp,yp,x,y,sxy,dx,dy);
        end
        intsum = intensity .* sdata;
        intsum = sum(intsum(:));
        
        sumsum = intensity.*intensity;
        sumsum = sum(sumsum(:)) ;
        
        sumx = (double(xp-0.5)).*intensity.*sdata;
        sumx = sum(sumx(:));
        
        sumy = (double(yp-0.5)).*intensity.*sdata;
        sumy = sum(sumy(:));
        
        if intsum <= 0 || sumsum == 0
            [x0(1,it),y0(1,it),N0(1,it)] = deal(-1,-1,-1);
        else
            x0(1,it) = double(sumx) / double(intsum);
            y0(1,it) = double(sumy) / double(intsum);
            N0(1,it) = double(intsum) / double(sumsum); 
            if location_out_of_ROI([x0(1,it),y0(1,it)],...
                    [xp_min,xp_max,yp_min,yp_max],[dx,dy])
                       [x0(1,it),y0(1,it),N0(1,it)] = deal(-1,-1,-1);
            end      
        end
        dist(1,it) = sqrt( (x-x0(1,it))^2 + (y-y0(1,it))^2 );
        x = x0(1,it); y = y0(1,it); 
        tmp = dist(1,it);
        if x0(1,it) == -1
            tmp = tol-1;
        end
        it = it+1; 
    end      
end

if it >1
    x0 = x0(1,it-1); 
    y0 = y0(1,it-1);
    if ndims(data) == 3
        z0 = z0(1,it-1);
    end
    N0 = N0(1,it-1);
    dist = dist(1,it-1);
    
end

%% error computation
if ndims(data) == 3
    err0 = N0*double(intensity) ...
        - double(data(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max));
    err0 = err0.^2;
    err0 = sqrt(sum(err0(:)));    
elseif ndims(data) == 2
    err0 = N0*double(intensity) - double(data(xp_min:xp_max,yp_min:yp_max));
    err0 = err0.^2;
    err0 = sqrt(sum(err0(:)));
end

%% convert PSF prefactor into integrated intensity
if ndims(data) == 3
    switch params.psfType
        case 'integratedGaussian'
           N0 = 8 * N0; 
        case 'integratedGaussianStdZ'
           x = (floor(x0 - 3*sxy): ceil(x0 + 3*sxy )  );
           y = (floor(y0 - 3*sxy): ceil(y0 + 3*sxy ) );
           z = (floor(z0 - 3*sz): ceil(z0 + 3*sz )  );
           [yy,xx,zz] = meshgrid(y,x,z);
           xx = ceil(xx)-0.5; 
           yy = ceil(yy)-0.5;
           zz = round(zz);
           Itot = intensity_integrated_gaussian3D_stdZ(...
               xx,yy,zz,x0,y0,z0,sxy,sz,1,1,1);
           N0 = sum(Itot(:))*N0; 
        case 'gaussian'
           x = (floor(x0 - 3*sxy): ceil(x0 + 3*sxy )  );
           y = (floor(y0 - 3*sxy): ceil(y0 + 3*sxy )  );
           z = (floor(z0 - 3*sz): ceil(z0 + 3*sz )  );
           [yy,xx,zz] = meshgrid(y,x,z);
           xx = ceil(xx)-0.5; 
           yy = ceil(yy)-0.5;
           zz = round(zz);
           Itot = intensity_gaussian3D(xx,yy,zz,x0,y0,z0,sxy,sz,1,1,1);
           N0 = sum(Itot(:))*N0;       
    end
else
    switch params.psfType
        case 'integratedGaussian'
           N0 = 4 * N0; 
        case 'integratedGaussianStdZ'
           N0 = 4 * N0;  
        case 'gaussian'
           x = (floor(x0 - 3*sxy): ceil(x0 + 3*sxy )  );
           y = ( floor(y0 - 3*sxy): ceil(y0 + 3*sxy )  );
           [yy,xx] = meshgrid(y,x);
           xx = ceil(xx)-0.5; 
           yy = ceil(yy)-0.5;
           Itot = intensity_gaussian2D(xx,yy,x0,y0,sxy,1,1);
           N0 = sum(Itot(:))*N0;    
    end
end

end

function res = location_out_of_ROI(location,boundaries,vox_size)
    %test if the location found by the mask is within the ROI
    %boundaries in pizel units
    %x0 y0 z0 refer to positions with origin on the edge of the pixel.
 
    if numel(location) == 3
        x0 = location(1);
        y0 = location(2);
        z0 = location(3);
        xp_min = boundaries(1);
        xp_max = boundaries(2);
        yp_min = boundaries(3);
        yp_max = boundaries(4);
        zp_min = boundaries(5);
        zp_max = boundaries(6);
        dx = vox_size(1);
        dy = vox_size(2);
        dz = vox_size(3);
    else
        x0 = location(1);
        y0 = location(2);
        xp_min = boundaries(1);
        xp_max = boundaries(2);
        yp_min = boundaries(3);
        yp_max = boundaries(4);
        dx = vox_size(1);
        dy = vox_size(2);
    end
    
    res = 0;    
    if x0/dx <xp_min-1
        res = 1;
        return
    end
    if x0/dx >xp_max
        res = 1;
        return
    end
    
    if y0/dy <yp_min -1
        res = 1;
        return
    end
    if y0/dy >yp_max
        res = 1;
        return
    end
    
    if numel(location) == 3
        if z0/dz <zp_min -1
            res = 1;
            return
        end
        if z0/dz >zp_max
            res = 1;
            return
        end
    end
    
end