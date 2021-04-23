function [new_ctr,ROIlimits] = compute_ROI_boundaries2(stackSize,...
    spot_ctr,cutwidth,thickness,ROIsize)

%check input consistency
if numel(stackSize) == 3 && numel(spot_ctr) < 3
    disp('cannot compute ROI boundaries, incorrect ROI center coordinates');
end

%compute half length of ROI square/cube
if strcmp(ROIsize,'small')
    halflength = cutwidth(1);
    if numel(stackSize) == 3
        halflength(2) = cutwidth(2);
    end
elseif strcmp(ROIsize,'large')
    halflength = cutwidth(1)+thickness;
    if numel(stackSize) == 3
        halflength(2) = cutwidth(2)+thickness;
    end
end

%input size
if numel(stackSize) == 3
    nx = stackSize(1);
    ny = stackSize(2);
    nz = stackSize(3);
    xc = spot_ctr(1); yc = spot_ctr(2); zc = spot_ctr(3);
else
    nx = stackSize(1);
    ny = stackSize(2);
    xc = spot_ctr(1); yc = spot_ctr(2); 
end

%compute ROI limits
xmin = ceil(xc - halflength(1));
xmin = max(1,xmin);
ymin = ceil(yc - halflength(1));
ymin = max(1,ymin);

xmax = ceil(xc + halflength(1));
xmax = min(nx,xmax);
ymax = ceil(yc + halflength(1));
ymax = min(ny,ymax);

ROIlimits = [xmin,ymin;xmax,ymax];

if numel(stackSize) == 3
    zmin = round(zc - halflength(2));
    zmin = max(1,zmin);
    
    zmax = round(zc + halflength(2));
    zmax = min(nz,zmax);
    
    ROIlimits = [ROIlimits,[zmin;zmax]];
end

%compute coordinates of spot center in new region
new_ctr(1) = spot_ctr(1) - xmin + 1;
new_ctr(2) = spot_ctr(2) - ymin + 1;
if numel(stackSize) == 3
    new_ctr(3) = spot_ctr(3) - zmin + 1;
end

end