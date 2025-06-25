% adaptive script
% given an image with spots and and an image of the same size of masks, it
% computes thresholds local to each ROI in the mask and detects spots in
% each mask. 
% (Spots outside of masks are tossed).

% MAKE SURE THAT THE FOLLOWING FOLDERS ARE ON THE PATH:
% AIRLOCALIZE/subfunctions
% AIRLOCALIZE/subfunctions/iniconfig

%% files paths for:
% - the data to be localized (imgName)
% - the masks used to generate locally relevant tresholds (maskName)
imgName = '/Users/lionnt01/Documents/junk/3D_stack_example.tif';
maskName = '/Users/lionnt01/Documents/junk/3D_stack_example_mask_lbl.tif';


imgName = '/Users/lionnt01/Documents/data/Nestor_for_Tim/AIRLOCALIZE_test/z_stacks_channel1/02112025_Cen-ZF_select_HCT116_NT-GFP_Cen7-ZF1g_1_MMStack_Pos2.ome.tif - C=0-1.tif';
maskName = '/Users/lionnt01/Documents/data/Nestor_for_Tim/AIRLOCALIZE_test/nuclei_masks/C1-02112025_Cen-ZF_select_HCT116_NT-GFP_Cen7-ZF1g_1_MMStack_Pos1_MIP_cp_masks.tif';
%% settings
% key parameter: size of the psf in pixel units 
% should be scalar for 2D image data / should be 2D for 3D z-stacks ([sigma_xy, sigma_z])
psfSigma = [1.3,2];
psfSigma = [2.3,3];

% key parameter: factor that will be multiplied to the std of the mask
% (recommended: 6)
threshFactor = 12;

% options
% if = 1, eliminates spots in the background region; set to 0 to keep (default: 1)
eliminateBackgroundSpots = 1;

% pixel value that marks the backgorund spots. If set to NaN, the region
% with minimal pixel value will be eliminated (default: 0)
backgroundID = 0;

% size of the pixel padding around each mask in the cropped image (default:
% 0)
paddingSize = 0;

% setup the data and parameters objects
alData1m = airLocalizeData();
params1m = airLocalizeParams();
params1m.reset();
params1m.threshUnits = 'adaptive';
params1m.threshLevel = 12;
params1m.psfSigma = psfSigma;

img = timtiffread(imgName);
mask = timtiffread(maskName);
%%
start0 = tic;
smooth = smooth_image_and_subtract_background9(img,params1m,'mask',mask);
t = toc(start0)
%% 
start2 = tic;
smooth = smooth_image_and_subtract_background9(img,params1m);
t = toc(start2)

%%
mm = unique(mask);
start1 = tic;
for i=1:numel(mm)
    %x = find(mask == mm(i));
    [croppedImg, dx, dy, dz] = cropImageBasedOnMask(img, mask, mm(i), 0);
    croppedMask = cropImageBasedOnMask(mask, mask, mm(i), 0);
end
toc(start1)

function [croppedImg, dx, dy, dz] = cropImageBasedOnMask(img, mask, mID, paddingSize)
% CROP_MASK Crops a region of 'img' based on a labeled region in 'mask'.
%
% [croppedImg, dx, dy, dz] = crop_mask(img, mask, mID, paddingSize)
%
% - img: 2D or 3D image (e.g., [nx, ny] or [nx, ny, nz])
% - mask: same size as img or a 2D mask if img is 3D
% - mID: integer value in mask to crop around
% - paddingSize: number of voxels to pad in each direction
%
% - croppedImg: cropped subregion of img
% - dx, dy, dz: offset indices into the original image

    % Validate inputs
    if nargin < 4
        paddingSize = 0;
    end

    nd = ndims(img);

    is3D = nd == 3;
    if size(img,1) ~= size(mask,1) || size(img,2) ~= size(mask,2)
        error('Mask must be the same size as image');
    end
    if is3D
        if ~ismatrix(mask)
            if size(img,3) ~= size(mask,3) 
                error('Mask must be 2D or the same size as img');
            end
        end
    end

    % Logical mask of the selected region
    region = (mask == mID);

    % Find bounding box
    if is3D
        if ismatrix(mask)
            % Convert ROI linear indices to subscript indices (i,j)
            [x, y] = ind2sub(size(region), find(region));
            z = [1,size(img,3)]; % dummy array of z values that encompasses the whole image
        else
            [x, y, z] = ind2sub(size(img), find(region));
        end
    else
        [x, y] = ind2sub(size(region), find(region));
    end
    xmin = max(min(x) - paddingSize, 1);
    xmax = min(max(x) + paddingSize, size(img, 1));
    ymin = max(min(y) - paddingSize, 1);
    ymax = min(max(y) + paddingSize, size(img, 2));

    if is3D
        zmin = max(min(z) - paddingSize, 1);
        zmax = min(max(z) + paddingSize, size(img, 3));
        croppedImg = img(xmin:xmax, ymin:ymax, zmin:zmax);
        dz = zmin;
    else
        croppedImg = img(xmin:xmax, ymin:ymax);
        dz = []; % no z offset for 2D
    end

    % Return offsets
    dx = xmin;
    dy = ymin;
end

