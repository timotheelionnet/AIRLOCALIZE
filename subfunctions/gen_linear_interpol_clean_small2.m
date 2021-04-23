function [pl,stack_bg_corr,new_ctr,ROIlimits,plmean,plfit,bg] = ...
gen_linear_interpol_clean_small2(alData,spotCenter,cutwidth,thickness,ROIsize)

%from a 3D (2D) stack 'stack' that contains a spot centered around the xc,yc,zc
%position (in pixels)
%returns another stack pl  
%values of pl are equal to the planar (3D linear) interpolation of
%the background surrounding the ROI (defined as a square of cutwidth pixels around the spot). 
%The background interpolation is done in a region of 'thickness' pixels from the boundary of the ROI
%ROIsize is an option: 'small' (default) is for the integrated function
%'large' generates a larger plane for debugging purpose.

%selects the points that are going to be used for the background correction
%i.e. located within thickness pixels of the ROI of the spot
%(ROI = square of size 2*cutwidth+1 size centered on spot_ctr =
%[xc,yc,zc]([xc,yc] for 2D)
%and stores them as a 4 columns array, resp [x y z I] ([x y I] for 2D
%all variables in pixel units
%ROIsize = 'large': output images include the extra layer used for background
%calculation; used for debugging
%ROIsize = 'small': output images exclude the extra layer used for background calculation
%version used for the integrated function

%outputs: 
%pl = the linear interpolation cropped over the ROI;
%stack_bg_corr = the original image corrected for background, cropped over the ROI   
%new_Ctr = [xc2,yc2,zc2] or [xc2,yc2] (2D) : the coordinates xc yc and zc in the generated images 
%ROIlimits = [xmin,ymin,zmin] or [xmin,ymin] (2D) = the coordinates of the ROI corner in the original image 
%ROIsize = 'large': ROI includes the extra layer used for background
%calculation;
%ROIsize = 'small': ROI excludes the extra layer used for background calculation
%plmean = the average background level over the ROI
%plfit: the fit parameters plfit = [a,b,c, d] where Ibackground = ax + by + cz + d
%plfit = [a,b,c] where Ibackground = ax + by + c for 2D
%a pixel [i,j,k], [x,y,z] = [i-0.5,j-0.5,k]; origin relative to original
%image
if ~alData.isMovie
    numDim = ndims(alData.img);
else
    numDim = ndims(alData.img) - 1;
    if numDim < 2
        disp('image dimensionality is not compatible with movie mode');
        return
    end
end
cutwidth_xy = cutwidth(1); 
xc = spotCenter(1);
yc = spotCenter(2);

if numDim==3 
    if numel(cutwidth)<2
        disp(['3D stack selected but cutwidth parameter has only 1 element;'...
            ' ensure that fit entry is compatible with 3D data']);
        [pl,stack_bg_corr,new_ctr,ROIlimits,plmean,plfit,bg] ...
            = deal(repmat([],1,7));
        return
    end
    cutwidth_z = cutwidth(2);
    if numel(spotCenter)<3
        disp(['3D stack selected but spot center has only 2 coordinates;'...
            ' ensure that numdim is set to 3.']);
        [pl,stack_bg_corr,new_ctr,ROIlimits,plmean,plfit,bg] ...
            = deal(repmat([],1,7));
        return
    end
    zc = spotCenter(3);
    bg = generate_bg_region_3D2(...
        alData,xc,yc,zc,cutwidth_xy,cutwidth_z,thickness);
else
    bg = generate_bg_region_2D2(alData,xc,yc,cutwidth_xy,thickness);
end

if isempty(bg)
    if numDim==3 
        plfit = zeros(1,4);
        new_ctr = [xc,yc,zc];
        ROIlimits = [1,1,1];
    else
        plfit = zeros(1,3);
        new_ctr = [xc,yc,0];
        ROIlimits = [1,1];
    end
    pl = [];
    stack_bg_corr = [];
    plmean = 0;
    return    
end

%fits the intensity distribution in these points to a 4D-plane
%equation ax+by+cz+d = intensity -- plfit = [a,b,c,d];
%all variables in pixel units
if numDim==3 
    plfit = fit_to_4D_plane2(bg);
else
    plfit = fit_to_3D_plane2(bg);
end

%generate an image with same size as the original,
%where the intensity values in the ROI are equal to the result of the
%plane fitting.
%the rest of the image is zeros.
%ROIsize can be set to 'large' for debugging purpose
%all variables in pixel units
if numDim ==3
    [pl,stack_bg_corr,new_ctr,ROIlimits] = ...
        final_4D_plane_small2(alData,xc,yc,zc,cutwidth_xy,cutwidth_z,...
        thickness,plfit,ROIsize);
else
    [pl,stack_bg_corr,new_ctr,ROIlimits] = ...
        final_3D_plane_small2(alData,xc,yc,cutwidth_xy,...
        thickness,plfit,ROIsize);
end
plmean = mean(pl(:));

end
    
function bg = generate_bg_region_3D2(alData,xcenter,ycenter,zcenter,...
    cutwidth_xy,cutwidth_z,thickness)

%%%%%%all variables in pixel units!
%from an image img with a spot located at xc,yc,zc 
%defines a region around the spot in order to compute the background intensity.
%It selects points located within 'thickness' pixels from the ROI.
%(ROI = cube of size 2*cutwidth+1 centered on xc,yc,zc).
%It returns bg: the array of points formatted as 4 cols: (x,y,z,intensity)
%the returned x,y,z values refer to the position of the center of the pixel (in pixel units),
%with the corner of the image defined as 0,0,0.

%%
thickness = floor(abs(thickness)); 
cutwidth_xy = floor(abs(cutwidth_xy));
cutwidth_z = floor(abs(cutwidth_z));

[nx,ny,nz] = size(alData.img);

npts = 4*thickness*( 2*cutwidth_xy + thickness + 1 )*( 2*cutwidth_z + 1 ) ...
+ 2*thickness*( 2*(thickness+cutwidth_xy) + 1 )*( 2*(thickness+cutwidth_xy) + 1 );

bg = zeros(npts+1,4);

%pixel # of the center of the region 
xc = ceil(xcenter); 
xc = max(1,min(xc,nx));
yc = ceil(ycenter); 
yc = max(1,min(yc,ny));
zc = ceil(zcenter);
zc = max(1,min(zc,nz));


%%
x2 = xc - cutwidth_xy;
x2 = max(x2,1);
y2 = yc - cutwidth_xy;
y2 = max(y2,1);
z2 = zc - cutwidth_z;
z2 = max(z2,1);

x3 = xc + cutwidth_xy;
x3 = min(x3,nx);
y3 = yc + cutwidth_xy;
y3 = min(y3,ny);
z3 = zc + cutwidth_z;
z3 = min(z3,nz);

if x2==1
    x1 = x2+1;
else
    x1 = ceil(xc) - cutwidth_xy - thickness;
    x1 = max(x1,1);
end

if y2==1
    y1 = y2+1;
else
    y1 = ceil(yc) - cutwidth_xy - thickness;
    y1 = max(y1,1);
end

if z2==1
    z1 = z2+1;
else
    z1 = ceil(zc) - cutwidth_z - thickness;
    z1 = max(z1,1);
end


if x3==nx
    x4 = x3-1;
else
    x4 = ceil(xc) + cutwidth_xy + thickness;
    x4 = min(x4,nx);
end

if y3==ny
    y4 = y3-1;
else
    y4 = ceil(yc) + cutwidth_xy + thickness;
    y4 = min(y4,ny);
end

if z3==nz
    z4 = z3-1;
else
    z4 = ceil(zc) + cutwidth_z + thickness;
    z4 = min(z4,nz);
end


%%
k=1;
for zpix = z2:z3
    for xpix = x1:x4
        for ypix = y1:y2-1
            bg(k,1) = xpix; 
            bg(k,2) = ypix;
            bg(k,3) = zpix;
            bg(k,4) = alData.img(xpix,ypix,zpix);
            k = k+1;
        end
    end

    for xpix=x1:x2-1
        for ypix=y2:y3
            bg(k,1) = xpix;
            bg(k,2) = ypix;
            bg(k,3) = zpix;
            bg(k,4) = alData.img(xpix,ypix,zpix);
            k = k+1;
        end
    end

    for xpix=x3+1:x4
        for ypix=y2:y3
            bg(k,1) = xpix;
            bg(k,2) = ypix;
            bg(k,3) = zpix;
            bg(k,4) = alData.img(xpix,ypix,zpix);
            k=k+1;
        end
    end

    for xpix=x1:x4
        for ypix=y3+1:y4
            bg(k,1) = xpix;
            bg(k,2) = ypix;
            bg(k,3) = zpix;
            bg(k,4) = alData.img(xpix,ypix,zpix);
            k=k+1;
        end
    end
end

for zpix = z1:z2-1
    for xpix = x1:x4
        for ypix = y1:y4
            bg(k,1) = xpix;
            bg(k,2) = ypix;
            bg(k,3) = zpix;
            bg(k,4) = alData.img(xpix,ypix,zpix);
            k=k+1;
        end
    end
end

for zpix = z3+1:z4
    for xpix = x1:x4
        for ypix = y1:y4
            bg(k,1) = xpix;
            bg(k,2) = ypix;
            bg(k,3) = zpix;
            bg(k,4) = alData.img(xpix,ypix,zpix);
            k=k+1;
        end
    end
end

%everything is in pixel coordinates so far (i.e. origin @ [1,1]). 
%I now transfer this into space coordinates (origin @ [0,0]). This means
%that pix i,j is referred at by its center located at [i-0.5,j-0.5]
bg = bg(1:k-1,:);
bg(:,1) = bg(:,1)-0.5;
bg(:,2) = bg(:,2)-0.5;
bg(:,3) = bg(:,3);

end

function bg = generate_bg_region_2D2(alData,xcenter,ycenter,cutwidth_xy,thickness)

%%%%%%all variables in pixel units!
%from an image img with a spot located at xc,yc 
%defines a region around the spot in order to compute the background intensity.
%It selects points located within 'thickness' pixels from the ROI.
%(ROI = square of size 2*cutwidth+1 centered on xc,yc).
%It returns bg: the array of points formatted as 3 cols: (x,y,intensity)
%the returned x,y values refer to the position of the center of the pixel (in pixel units),
%with the corner of the image defined as 0,0.

%%
thickness = floor(abs(thickness)); 
cutwidth_xy = floor(abs(cutwidth_xy));
if ~alData.isMovie
    [nx,ny] = size(alData.img);
else
    [nx,ny] = size(alData.img(:,:,alData.curFrame));
end
npts =  2*thickness*( 2*(thickness+cutwidth_xy) + 1 );
bg = zeros(npts+1,3);

%pixel # of the center of the region 
xc = ceil(xcenter); 
xc = max(1,min(xc,nx));
yc = ceil(ycenter); 
yc = max(1,min(yc,ny));

%%
x2 = xc - cutwidth_xy;
x2 = max(x2,1);
y2 = yc - cutwidth_xy;
y2 = max(y2,1);

x3 = xc + cutwidth_xy;
x3 = min(x3,nx);
y3 = yc + cutwidth_xy;
y3 = min(y3,ny);

x1 = ceil(xc) - cutwidth_xy - thickness;
x1 = max(x1,1);
y1 = ceil(yc) - cutwidth_xy - thickness;
y1 = max(y1,1);

x4 = ceil(xc) + cutwidth_xy + thickness;
x4 = min(x4,nx);
y4 = ceil(yc) + cutwidth_xy + thickness;
y4 = min(y4,ny);


%%
k=1;
for xpix = x1:x4
    for ypix = y1:y2-1
        bg(k,1) = xpix; 
        bg(k,2) = ypix;
        if ~alData.isMovie
            bg(k,3) = alData.img(xpix,ypix);
        else
            bg(k,3) = alData.img(xpix,ypix,alData.curFrame);
        end
        k = k+1;
    end
end

for xpix=x1:x2-1
    for ypix=y2:y3
        bg(k,1) = xpix;
        bg(k,2) = ypix;
        if ~alData.isMovie
            bg(k,3) = alData.img(xpix,ypix);
        else
            bg(k,3) = alData.img(xpix,ypix,alData.curFrame);
        end
        k = k+1;
    end
end

for xpix=x3+1:x4
    for ypix=y2:y3
        bg(k,1) = xpix;
        bg(k,2) = ypix;
        if ~alData.isMovie
            bg(k,3) = alData.img(xpix,ypix);
        else
            bg(k,3) = alData.img(xpix,ypix,alData.curFrame);
        end
        k=k+1;
    end
end

for xpix=x1:x4
    for ypix=y3+1:y4
        bg(k,1) = xpix;
        bg(k,2) = ypix;
        if ~alData.isMovie
            bg(k,3) = alData.img(xpix,ypix);
        else
            bg(k,3) = alData.img(xpix,ypix,alData.curFrame);
        end
        k=k+1;
    end
end


%everything is in pixel coordinates so far (i.e. origin @ [1,1]). 
%I now transfer this into space coordinates (origin @ [0,0]). This means
%that pix i,j is referred at by its center located at [i-0.5,j-0.5]
bg = bg(1:k-1,:);
bg(:,1) = bg(:,1)-0.5;
bg(:,2) = bg(:,2)-0.5;

end

function x = fit_to_3D_plane2(data)
%fits the intensity I vs. (x,y) to a 3D-hyperplane with equation I = ax + by + c.
%data should be formatted as a 3 columns: x,y,intensity.
%where x,y refer to the physical coordinates (in pix units) 
%of the center of the pixel (origin 0,0 at the corner of the image).
%returns a vector x = [a,b,c]



npts = size(data,1);
sx = sum(data(:,1));
sy = sum(data(:,2));
sI = sum(data(:,3));
sxx = sum( data(:,1).^2 );
syy = sum( data(:,2).^2 );
sxy = sum( data(:,1).*data(:,2) );
sIx = sum( data(:,3).*data(:,1) );
sIy = sum( data(:,3).*data(:,2) );

fitmat = [sxx sxy sx; sxy syy sy; sx sy npts];
DD = det(fitmat);

%if matrix inversion is impossible, I just return the average intensity as
%the result

if (DD==0) 
    x = [0;0;mean(data(:,3))];
else  
    v = [sIx;sIy;sI];
    x = fitmat\v;
end

clear('data',...
    'sx','sy','sI',...
    'sxx','syy',...
    'sxy','sxz',...
    'sIx','sIy',...
    'fitmat','DD','v','npts');
end

function x = fit_to_4D_plane2(data)
%fits the intensity I vs. (x,y,z) to a 4D-hyperplane with equation I = ax + by + cz +d.
%data should be formatted as a 4 columns: x,y,z,intensity.
%where x,y,z refer to the physical coordinates (in pix units) 
%of the center of the pixel (origin 0,0 at the corner of the image).
%returns a vector x = [a,b,c,d]
data = double(data);

npts = size(data,1);
sx = sum(data(:,1));
sy = sum(data(:,2));
sz = sum(data(:,3));
sI = sum(data(:,4));
sxx = sum( data(:,1).*data(:,1) );
syy = sum( data(:,2).*data(:,2) );
szz = sum( data(:,3).*data(:,3) );
sxy = sum( data(:,1).*data(:,2) );
sxz = sum( data(:,1).*data(:,3) );
syz = sum( data(:,2).*data(:,3) );
sIx = sum( data(:,4).*data(:,1) );
sIy = sum( data(:,4).*data(:,2) );
sIz = sum( data(:,4).*data(:,3) );

fitmat = [sxx sxy sxz sx; sxy syy syz sy; sxz syz szz sz; sx sy sz npts];
DD = det(fitmat);
cond(fitmat);
%if matrix inversion is impossible, I just return the average intensity as
%the result
if (DD==0) 
    x = [0;0;0;mean(data(:,4))];
else
    v = [sIx;sIy;sIz;sI];
    x = fitmat\v;
end

clear('data',...
    'sx','sy','sz','sI',...
    'sxx','syy','szz',...
    'sxy','sxz','syz',...
    'sIx','sIy','sIz',...
    'fitmat','DD','v','npts');
end

function  [pl, stack_bg_corr,new_ctr,ROIlimits] = ...
    final_4D_plane_small2(alData,xc,yc,zc,cutwidth_xy,cutwidth_z,...
    thickness,plfit,ROIsize)
    %function that generates an image pl
    %the image has the size of an ROI defined as within 'cutwidth' pixels 
    %from the center of a spot xc,yc from the mother image.
    %in the ROI, the pixels have values equal to the planar interpolation of
    %the background surrounding the ROI (parameters of the plane in 'plfit'
    %vector).
    %the size of the ROI can be made larger (including the background around
    %it) when selecting the 'large' option.

    [new_ctr,ROIlimits] = compute_ROI_boundaries2(...
        size(alData.img),[xc,yc,zc],[cutwidth_xy,cutwidth_z],thickness,ROIsize);
    
    xmin = ROIlimits(1,1); 
    xmax = ROIlimits(2,1); 
    ymin = ROIlimits(1,2); 
    ymax = ROIlimits(2,2);
    zmin = ROIlimits(1,3); 
    zmax = ROIlimits(2,3);
    
    pl = zeros(xmax-xmin+1,ymax-ymin+1,zmax-zmin+1);
    stack_bg_corr = zeros(xmax-xmin+1,ymax-ymin+1,zmax-zmin+1);
    
    for xpix=xmin:xmax
        for ypix=ymin:ymax
            for zpix=zmin:zmax

            %4D plane parameters 'plfit' correspond to the physical coordinates (in pix units) 
            %of the center of the pixel (origin 0,0 at the corner of the image).
            pl(xpix-xmin+1,ypix-ymin+1,zpix-zmin+1) = ...
                plfit(1)*(xpix-0.5) + plfit(2)*(ypix-0.5) ...
                + plfit(3)*zpix + plfit(4);
            stack_bg_corr(xpix-xmin+1,ypix-ymin+1,zpix-zmin+1) = ...
                double(alData.img(xpix,ypix,zpix)) - ...
                pl(xpix-xmin+1,ypix-ymin+1,zpix-zmin+1);
            end
        end
    end
clear('stack','nx','ny','nz',...
    'xpix','ypix','zpix');
    
end

function  [pl, stack_bg_corr,new_ctr,ROIlimits] = ...
    final_3D_plane_small2(alData,xc,yc,cutwidth_xy,...
    thickness,plfit,ROIsize)
    %function that generates an image pl
    %the image has the size of an ROI defined as within 'cutwidth' pixels 
    %from the center of a spot xc,yc from the mother image.
    %in the ROI, the pixels have values equal to the planar interpolation of
    %the background surrounding the ROI (parameters of the plane in 'plfit'
    %vector).
    %the size of the ROI can be made larger (including the background around
    %it) when selecting the 'large' option.
    if ~alData.isMovie
        imSize = size(alData.img);
    else
        imSize = size(alData.img(:,:,alData.curFrame));
    end
    
    [new_ctr,ROIlimits] = compute_ROI_boundaries2(...
        imSize,[xc,yc],cutwidth_xy,thickness,ROIsize);
    
    xmin = ROIlimits(1,1); 
    xmax = ROIlimits(2,1); 
    ymin = ROIlimits(1,2); 
    ymax = ROIlimits(2,2);    
    
    pl = zeros(xmax-xmin+1,ymax-ymin+1);
    stack_bg_corr = zeros(xmax-xmin+1,ymax-ymin+1);
    
    
    for xpix=xmin:xmax
        for ypix=ymin:ymax
            %4D plane parameters 'plfit' correspond to the physical coordinates (in pix units) 
            %of the center of the pixel (origin 0,0 at the corner of the image).
            pl(xpix-xmin+1,ypix-ymin+1) = ...
                plfit(1)*(xpix-0.5) + plfit(2)*(ypix-0.5) + plfit(3);
            if ~alData.isMovie
                stack_bg_corr(xpix-xmin+1,ypix-ymin+1) = ...
                    double(alData.img(xpix,ypix)) ...
                    - pl(xpix-xmin+1,ypix-ymin+1);
            else
                stack_bg_corr(xpix-xmin+1,ypix-ymin+1) = ...
                    double(alData.img(xpix,ypix,alData.curFrame)) ...
                    - pl(xpix-xmin+1,ypix-ymin+1);
            end
        end
    end
    
clear('stack','nx','ny','nz',...
    'xpix','ypix','zpix');
    
end



    
    
    
    