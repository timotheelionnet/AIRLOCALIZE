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

% nestor
imgName = ['/Users/lionnt01/Documents/data/Nestor_for_Tim/AIRLOCALIZE_test/z_stacks_channel1/',...
    '02112025_Cen-ZF_select_HCT116_NT-GFP_Cen7-ZF1g_1_MMStack_Pos1.ome.tif - C=0-1.tif'];
maskName = ['/Users/lionnt01/Documents/data/Nestor_for_Tim/AIRLOCALIZE_test/nuclei_masks/',...
    'C1-02112025_Cen-ZF_select_HCT116_NT-GFP_Cen7-ZF1g_1_MMStack_Pos1_MIP_cp_masks.tif'];

% minghan
imgName = '/Users/lionnt01/Documents/data/v1_adaptive/img_channel001_position002_time000000000.tif';
maskName = ['/Users/lionnt01/Documents/data/v1_adaptive_mask/',...
    'MAX_img_channel001_position002_time000000000_cp_masks.tif'];
%% settings
% key parameter: size of the psf in pixel units 
% should be scalar for 2D image data / should be 2D for 3D z-stacks ([sigma_xy, sigma_z])
psfSigma = [1.3,2];
%psfSigma = [2.3,3]; % nestor setting

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
%%
img = timtiffread(imgName);
mask = timtiffread(maskName);
%%
loList = 0.2:0.1:3;
hiList = 0.2:0.1:3;
numRegions = zeros(numel(loList),numel(hiList));
for i = 1:numel(loList)
    disp(['i = ',num2str(i),'/',num2str(numel(loList)),'...']);
    for j = 1:numel(hiList)
        params1m.filterLo = loList(i);
        params1m.filterHi = loList(j);
        smooth = 1000*smooth_image_and_subtract_background10(img,params1m,...
            'mask',mask,'paddingSize',5);
        t = toc(start0);
        s = smooth(:,:,11);
        m = max(s(mask == 14));
        s(s<=m | mask ~= 17) = 0;
        cc = bwconncomp(s);            % Find connected components
        numRegions(i,j) = cc.NumObjects; 
    end
end
save_as_tiff(smooth,'/Users/lionnt01/Documents/junk/testSmooth.tif');
%%
[l,m] = ind2sub(size(numRegions), find(numRegions == max(numRegions(:))));
loList(l)
hiList(m)

loList(13)
hiList(6)

%% 
start2 = tic;
smooth = smooth_image_and_subtract_background10(img,params1m);
t = toc(start2)
%%
start2 = tic;
smooth = smoothInMasks(img,mask,params1m.filterHi, params1m.filterLo, params1m.psfSigma(1),...
    1,0,5);
save_as_tiff(smooth,'/Users/lionnt01/Documents/junk/testSmooth.tif');
t = toc(start2)
% %%
% function smooth = smoothInMasks(img,mask,filterHi, filterLo, psfSigma_xy,...
%     eliminateBackgroundROIs,backgroundID,paddingSize)
% 
%     % replicate the mask in the third dim if it is two-dimensional and the
%     % image is a stack.
%     nd = ndims(img);
% 
%     % collect list of masks IDs
%     mList = collectMaskIds(mask,eliminateBackgroundROIs,backgroundID);
%     nm = numel(mList);
%     if ismember(0,mList)
%         disp('0 is part of the mask list');
%     end
% 
%     %% loop through masks
%     smooth = zeros(size(img));
%     t = zeros(nm,10);
%     %for i=1:nm 
%     for i=14:14     
%         % crop rectangle ROI around mask 
%         t1 = tic;
%         [croppedImg,dx,dy,dz,idxCrop] = cropImageBasedOnMask(...
%             img, mask, mList(i), paddingSize);
%         t(i,1) = toc(t1);
%         if ndims(mask) == ndims(img)
%             croppedMask = false(size(croppedImg));
%             croppedMask(idxCrop) = true;
%         else
%             % keep cropped mask 2D when possible for computation efficiency
%             croppedMask = mask(dx:dx+size(croppedImg,1)-1,...
%                 dy:dy+size(croppedImg,2)-1) == mList(i);
%         end
% 
%         % set the background around the mask to the value of the closest pixel
%         % within the image - this mitigates boundary artefacts.
%         t2 = tic;
%         croppedImg = padWithObjectValues(croppedImg, croppedMask, true);
%         %croppedImg = padWithObjectValuesLinIdx(croppedImg, idxCrop);
%         t(i,2) = toc(t2);
% 
%         % smooth the cropped image
%         t3 = tic;
%         cs = smoothImg(croppedImg, filterHi, filterLo, psfSigma_xy);
%         t(i,3) = toc(t3);
% 
%         % get mean & std of intensity over the mask region
%         t4 = tic;
%         %[~,s,m2,s2] = getMeanStdInMask(cs,croppedMask,true);
%         median(cs)
%         [mmm,s,m2,s2] = getMeanStdInMaskLinIdx(cs,idxCrop,'downSample',100)
%         t(i,4) = toc(t4);
% 
%         % subtract the mode within the mask, and divide by the estimate of
%         % the STD
%         croppedSmooth = zeros(size(cs));
%         croppedSmooth(idxCrop) = cs(idxCrop);
%         %croppedSmooth(idxCrop) = ( croppedSmooth(idxCrop) - max(0,m2) ) / ((s2+s)/2);
%         croppedSmooth(idxCrop) = ( croppedSmooth(idxCrop) - m2 ) / ((s2+s)/2);
%         croppedSmooth(isinf(croppedSmooth) | croppedSmooth<0) = 0; % avoid any infinites, negatives
% 
%         % reassign pixels in the large image to the smoothed/normalized
%         % values
%         t5 = tic;
%         if nd==2
%             [nx, ny] = size(cs);
%             smooth(dx:dx+nx-1,dy:dy+ny-1) = smooth(dx:dx+nx-1,dy:dy+ny-1) ...
%                 + croppedSmooth;
%         else
%             [nx, ny, nz] = size(cs);
%             smooth(dx:dx+nx-1,dy:dy+ny-1,dz:dz+nz-1) = ...
%                 smooth(dx:dx+nx-1,dy:dy+ny-1,dz:dz+nz-1) ...
%                 + croppedSmooth;
%         end
%         t(i,5) = toc(t5);
%     end
%     sum(t)
% end
% 
% %%
% function smooth = smoothImg(img, filterHi, filterLo, psfSigma_xy)
%     if ismember(ndims(img),[2,3])
% 
%         %factors 1.5 is a heuristic attempt to reproduce results
%         %from the Fourier filter given the same parameters
%         smooth = 1.5 * DOGfilter2(img, filterHi, filterLo * psfSigma_xy, [],[]);
%     else
%         disp('data type error: smoothing only possible on 2d images or 3d stacks');
%         smooth = 0;
%         return;
%     end
%     % if ismember(ndims(img),[2,3])
%     % 
%     %     %factors 1.5 is a heuristic attempt to reproduce results
%     %     %from the Fourier filter given the same parameters
%     %     smooth = 1.5 * DOGfilter3(img, psfSigma_xy, 3, [],[]);
%     % else
%     %     disp('data type error: smoothing only possible on 2d images or 3d stacks');
%     %     smooth = 0;
%     %     return;
%     % end
% end
% 
% %%
% function mList = collectMaskIds(mask,eliminateBackgroundROIs,backgroundID)
%     mList = unique(mask(:));
% 
%     % eliminate background spots if needed
%     if eliminateBackgroundROIs
%         if isnan(backgroundID)
%             mList = setdiff(mList,min(mList(:)));
%         else
%             mList = setdiff(mList,backgroundID);
%         end
%     end
% end
% 
% %%
% function [croppedImg, dx, dy, dz, idxCrop, idxImg] = cropImageBasedOnMask(img, mask, mID, paddingSize)
% % CROP_MASK Crops a region of 'img' based on a labeled region in 'mask'.
% %
% % [croppedImg, dx, dy, dz] = crop_mask(img, mask, mID, paddingSize)
% %
% % INPUT
% % - img: 2D or 3D image (e.g., [nx, ny] or [nx, ny, nz])
% % - mask: same size as img or a 2D mask if img is 3D
% % - mID: integer value in mask to crop around
% % - paddingSize: number of voxels to pad in each direction
% %
% % - croppedImg: cropped subregion of img
% % - dx, dy, dz: offset indices into the original image
% % idxCrop: list of linear indices of the voxels that map to the object relative to the
%     % cropped image
% % idxImg: list of linear indices of the voxels that map to the object
%     % relative the original image
% 
%     % Validate inputs
%     if nargin < 4
%         paddingSize = 0;
%     end
% 
%     % this should have been checked upstream so commenting it out
%     % % check size/dimensions compatibility between img and mask
%     % [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask);
%     % if errorFlag ==1
%     %     croppedImg = []; dx=[]; dy=[]; dz=[]; idxCrop=[]; idxImg = [];
%     %     return
%     % end
% 
%     nd = ndims(img);
%     is3D = nd == 3;
% 
%     % Logical mask of the selected region
%     idxImg = find(mask == mID);
% 
%     % Find bounding box
%     if is3D
%         [nx, ny, nz] = size(img);
%         if ismatrix(mask)
%             % Convert ROI linear indices to subscript indices (i,j)
%             [x, y] = ind2sub([nx, ny], idxImg);
%             z = [1,size(img,3)]; % dummy array of z values that encompasses the whole image
%         else
%             [x, y, z] = ind2sub([nx, ny, nz], idxImg);
%         end
%     else
%         [nx, ny] = size(img);
%         [x, y] = ind2sub([nx, ny], idxImg);
%     end
%     xmin = max(min(x) - paddingSize, 1);
%     xmax = min(max(x) + paddingSize, nx);
%     ymin = max(min(y) - paddingSize, 1);
%     ymax = min(max(y) + paddingSize, ny);
% 
%     if is3D
%         zmin = max(min(z) - paddingSize, 1);
%         zmax = min(max(z) + paddingSize, nz);
%         croppedImg = img(xmin:xmax, ymin:ymax, zmin:zmax);
%         dz = zmin;
%     else
%         croppedImg = img(xmin:xmax, ymin:ymax);
%         dz = []; % no z offset for 2D
%     end
% 
%     % Return offsets
%     dx = xmin;
%     dy = ymin;
% 
%     % compute linear indices relative to the whole image
%     if is3D && ismatrix(mask)
% 
%     end
% 
%     % compute linear indices relative to the mask
%     if is3D
%         if ismatrix(mask)
%             idxCrop = expand2DIndicesTo3D(...
%             [reshape(x - xmin + 1, numel(x),1), ...
%              reshape(y - ymin + 1, numel(y),1)],...
%             size(croppedImg,1),size(croppedImg,2),nz,1,'subscript');
%         else
%             idxCrop = sub2ind(size(croppedImg),...
%             x - xmin + 1,y - ymin + 1,z - zmin + 1);
%         end
%     else
%         idxCrop = sub2ind(size(croppedImg),...
%             x - xmin + 1,y - ymin + 1);
%     end
% 
%     % if is3D && ismatrix(mask)
%     %     idxImg = expand2DIndicesTo3D(idxImg,nx,ny,nz,1,'linear');
%     % end
%     % 
%     % % compute linear indices relative to the mask
%     % if is3D
%     %     if ismatrix(mask)
%     %         [x, y, z] = ind2sub([nx, ny, nz], idxImg);
%     %     end
%     %     x = x - xmin + 1;
%     %     y = y - ymin + 1;
%     %     z = z - zmin + 1;
%     %     idxCrop = sub2ind(size(croppedImg),x,y,z);
%     % else
%     %     x = x - xmin + 1;
%     %     y = y - ymin + 1;
%     %     idxCrop = sub2ind(size(croppedImg),x,y);
%     % end
% end
% 
% %%
% function idx3D = expand2DIndicesTo3D(idx2D,nx,ny,nz,squeezeResult,inputMode)
%     % from  an array of 2D linear indices (or x,y coordinates) relative to a matrix of size [nx,ny], 
%     % generates an array of linear indices of all the voxels of the 3D array
%     % (size [nx ny nz]) which 2D projection is a member of the list of 2D
%     % linear indices.
%     % if the input has size [ni nj], the output has size [ni nj nz]
%     % if the option squeezeResult has been selected, the dimensions of zie
%     % 1 will be squeezed to return the lowest possible dimension array.
% 
%     % INPUT
%     % idx2D either a list of linear indices (if inputMode = 'linear'), 
%         % or a 2D array or x y coordinates in 2 columns (if inputMode = 'subscript')
%     % nx,ny,nz is the size of the 3D array that the indices are replicated into
%     % squeezeResult is to ensure the dimensions of the result are squeezed.
% 
%     % inputMode can be either: 
%         % 'linear': idx2D has the format of linear indices
%         % or 
%         % 'subscript': idx2D is a 2 column matrix with x and y coordinates 
% 
%     % get list of x and y coordinates    
%     switch inputMode 
%         case 'linear'
%             % map idx2D array to a vector
%             s = size(idx2D);
%             idx2D  = idx2D(1:numel(idx2D));
% 
%             % Convert to subscript indices (i,j)
%             [i, j] = ind2sub([nx, ny], idx2D);
%         case 'subscript'
%             i = idx2D(:,1);
%             j = idx2D(:,2);
%             s = size(idx2D(:,1));
%     end
% 
%     % Expand to all k values
%     k = reshape(1:nz, 1, nz);  % [1 x nz]
%     i = reshape(i, [], 1);     % [N x 1]
%     j = reshape(j, [], 1);     % [N x 1]
% 
%     % Broadcast and get all (i,j,k)
%     I = repmat(i, 1, nz);  % [N x nz]
%     J = repmat(j, 1, nz);  % [N x nz]
%     K = repmat(k, length(i), 1);  % [N x nz]
% 
%     % Convert back to linear indices into [nx ny nz]
%     idx3D = sub2ind([nx ny nz], I(:), J(:), K(:));
% 
%     % recover original array format of the input
%     idx3D = reshape(idx3D,[s,nz]);
% 
%     % remove useless dimensions if option has been selected
%     if squeezeResult
%         idx3D = squeeze(idx3D);
%     end
% end
% 
% %%
% function paddedImg = padWithObjectValues(img, mask, maskVal)
% 
%     % this should have been checked upstream of this function call so
%     % commenting it out
%     % % check that size and dimensions of img and mask are compatible 
%     % [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask);
%     % if errorFlag == 1
%     %     paddedImg = [];
%     %     return
%     % end
% 
%     % Ensure mask is logical
%     mask = (mask == maskVal);
% 
%     % Compute distance transform and index of nearest object pixel
%     [~, idx] = bwdist(mask);
% 
%     % if mask is 2D and image is 3D, convert the 2D closest pixel in the
%     % mask into its 3D counterpart
%     if ndims(img) == 3 && ismatrix(mask)
%         squeezeResult = 0;
%         [nx,ny,nz] = size(img);
%         idx = expand2DIndicesTo3D(idx,nx,ny,nz,squeezeResult,'linear');
%     end
% 
%     % Initialize result
%     paddedImg = img;
% 
%     % Fill background pixels with values from nearest object pixels
%     paddedImg(~mask) = img(idx(~mask));
% end
% 
% 
% function [m,s,m2,s2,idxList] = getMeanStdInMaskLinIdx(img,idxList,varargin)
%     % computes summary statistics in an ROI of an img indexed by linIdx
% 
%     % INPUT
%     % img is the image
%     % idxList is the list of linear indices that define the ROI
% 
%     % OPTIONAL
%     % 'useMode': set to 1 to output the robust estimate of the mode for m2
%     % instead of the median (Default is zero, output m2 = median).
% 
%     % OUTPUT
%     % m: mean of pixels that belonw to the mask
%     % s: std of pixels that belonw to the mask
%     % m2: median of pixels that belonw to the mask
%     % s2: pseudo std of pixels that belonw to the mask
%         % computed as m2 - median( pixels with values below m2)
%     % idxList: list of linear indices of the image voxels that fall under
%         % the mask
% 
%     % the function handles img and mask having the same size, 
%     % or img being 3D and mask being 2D, in which case all img voxels that project
%     % vertically onto the ROI in the mask are included in the summary
%     % statistics.
% 
%     % gathe optional args
%     p = inputParser();
%     p.addParameter('useMode',0);
%     p.addParameter('downSample',1);
%     p.parse(varargin{:});
%     useMode = p.Results.useMode;
%     downSample = p.Results.downSample;
% 
%     % this should have been checked upstream of this function call so
%     % commenting it out
%     % % check size/dimensions compatibility between img and mask
%     % [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask);
%     % if errorFlag ==1
%     %     m = NaN; s = NaN;m2 = NaN; s2 = NaN;
%     %     return
%     % end
% 
%     % downsample the indices if needed to accelerate, convert to vector
%     img = img(idxList(1:ceil(downSample):numel(idxList))); 
%     m = mean(img);
%     s = std(img);
%     if useMode
%         % assing m2 to the mode if option is chosen
%         % number of bins used to compute the mode 
%         % (these bins will map the 20th 80th percentile interval)
%         nBins = 20;
%         m2 = estimateModeRobust(img,nBins);
%         s2 = m2 - median(img(img <= m2));
%     else
%         % assign m2 to median by default
%         m2 = prctile(img,[25,50])
%         s2 = m2(2) - m2(1);
%         m2 = m2(2);
%     end   
% end





