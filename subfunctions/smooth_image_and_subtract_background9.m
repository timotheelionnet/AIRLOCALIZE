function smooth = smooth_image_and_subtract_background9(img,params,varargin)
    % applies a bandpass filter to the image img
    % if image is a stack, each slice is 2D-filtered independently.
    
    % INPUT
    % img: the image or stack to be smoothed.
    % params: airlocalizeParams object
        % the function uses the following object properties
        % params.numDim: whether image is 2D or 3D
        % params.filterLo: low frequency cutoff in length in units of psf
            % Sigma
        % params.filterHi: hi frequency cutiff length in pixel units
        % params.psfSigma: lateral extension of the gaussian point spread function

    % OPTIONAL INPUT (name/value pairs):
    % 'mask': mask image - if populated, the smoothing is performed separately 
        % in each ROI of the mask. This allows adjusting for different
        % nuclei brighnesses for instance.
    % 'eliminateBackgroundROIs': if set to 1 (default), the ROI with 
        % ID equal to 'backgroundID' (i.e. the background) will be ignored; 
        % default backgroundID value is 0.
    % 'backgroundID': can be set to define the ID value of the background
        % ROI (default is 0). If 'backgroundID' is set to NaN, the background 
        % is picked as the ROI with minimum ID value in each mask.
    % 'paddingSize': the size in voxels surrounding each ROI during the
        % filtering step (default: 0).
    
    % OUTPUT
    % smooth: image or stack the same size as img which holds:
        % the smoothed data minus the average value of the smoothed image;
        % or if the mask option is selected:
        % the smoothed data minus the mask-wise pseudo-average value of the smoothed image, 
        % divided by the mask-wise pseudo-standard deviation of the
        % smoothed image. This is essentially a pseudo z-score.
        % The deviation from the canonical z-score formula provides
        % robustness against the (typical) non-normal distribution shapes
        % of fluorescence intensities.
        % smooth is set to zero in the pixels forming the background of the
        % mask.
    p = inputParser();
    addParameter(p,'mask',[]);
    addParameter(p,'eliminateBackgroundROIs',1);    
    addParameter(p,'backgroundID',0);    
    addParameter(p,'paddingSize',0);    
    
    parse(p,varargin{:});
    mask = p.Results.mask;
    eliminateBackgroundROIs = p.Results.eliminateBackgroundROIs;
    backgroundID = p.Results.backgroundID;
    paddingSize = p.Results.paddingSize;
    
    if isempty(mask)
        disp('smoothing image...'); 
        % classic bandpass smoothing filter applied to entire image
        smooth = smoothImg(img, params.numDim, params.filterHi, params.filterLo,...
            params.psfSigma(1));
    
        % subtract the background mean
        smooth = smooth - mean(smooth(:));
    else
        disp('smoothing image and normalizing intensity in masks...'); 
        smooth = smoothInMasks(img,mask,...
            params.filterHi, params.filterLo,params.psfSigma(1),...
            eliminateBackgroundROIs,backgroundID,paddingSize);
    end 
end

function smooth = smoothInMasks(img,mask,filterHi, filterLo, psfSigma_xy,...
    eliminateBackgroundROIs,backgroundID,paddingSize)

    % replicate the mask in the third dim if it is two-dimensional and the
    % image is a stack.
    nd = ndims(img);
    
    % collect list of masks IDs
    mList = collectMaskIds(mask,eliminateBackgroundROIs,backgroundID);
    nm = numel(mList);
    if ismember(0,mList)
        disp('0 is part of the mask list');
    end
    
    %% loop through masks
    smooth = zeros(size(img));
    for i=1:nm
        % crop rectangle ROI around mask 
        [croppedImg,dx,dy,dz] = cropImageBasedOnMask(...
            img, mask, mList(i), paddingSize);
        if nd == 3 
            if ismatrix(mask)
                % if image is 3D but mask is 2D, replicate the mask across
                % all the planes
                croppedMask = mask(dx:(dx+size(croppedImg,1)-1),...
                    dy:(dy+size(croppedImg,2)-1));
            else
                croppedMask = mask(dx:(dx+size(croppedImg,1)-1),...
                    dy:(dy+size(croppedImg,2)-1),...
                    dz:(dz+size(croppedImg,3)-1));
            end
        else
            croppedMask = mask(dx:(dx+size(croppedImg,1)-1),...
                    dy:(dy+size(croppedImg,2)-1));
        end
        
        % set the background around the mask to the value of the closest pixel
        % within the image - this mitigates boundary artefacts.
        croppedImg = padWithObjectValues(croppedImg, croppedMask, mList(i));
        
        % smooth the cropped image
        croppedSmooth = smoothImg(croppedImg, nd, filterHi, filterLo, psfSigma_xy);
        mList(i)
        min(croppedSmooth(:))
        max(croppedSmooth(:))
        % get mean & std of intensity over the mask region
        [~,s,m2,s2] = getMeanStdInMask(croppedSmooth,croppedMask,mList(i));
        
        % subtract the mode within the mask, and divide by the estimate of
        % the STD
        croppedSmooth = ( croppedSmooth - max(0,m2) ) / ((s2+s)/2);
        croppedSmooth(isinf(croppedSmooth)) = 0; % avoid any infinites
        croppedSmooth(croppedSmooth<0) = 0;
        
        % reassign pixels in the large image to the smoothed/normalized
        % values
        idxList = getIdxInImgMatchingMaskValue(...
                        croppedSmooth,croppedMask,mList(i));
        if nd==2
            [x,y] = ind2sub(size(croppedSmooth),idxList);
            x = x + dx -1; y = y + dy -1; 
            newIdxList = sub2ind(size(smooth),x,y);
        else
            [x,y,z] = ind2sub(size(croppedSmooth),idxList);
            x = x + dx -1; y = y + dy -1; z = z + dz -1;
            newIdxList = sub2ind(size(smooth),x,y,z);
        end
        smooth(newIdxList) = croppedSmooth(idxList);
    end
end

%%
function idx = getIdxInImgMatchingMaskValue(img,mask,maskVal)
% find indices in image img that map to the regions of the mask which have the
% pixel value maskVal. 
% Returns idx, a 1D array of the linear indices of the pixels. 
% Behavior is obvious if the sizes of img and mask are matching. 
% If img is 3D and mask is 2D, returns the indices of all the voxels of img 
% which x,y match maskVal when projected vertically onto mask.
    
    [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask);
    if errorFlag == 1
        idx = [];
        return
    end
    
    % trivial case when image and mask have the same sizes
    if isequal(size(img),size(mask))
        idx = find(mask==maskVal);
        return
    end

    % now is the case where image is 3D and mask is 2D
    idx2D = find(mask==maskVal);
    [nx, ny, nz] = size(img);
    squeezeResult = 1;
    idx = expand2DIndicesTo3D(idx2D,nx,ny,nz,squeezeResult); 
    idx = idx(1:numel(idx));
end

%%
function idx3D = expand2DIndicesTo3D(idx2D,nx,ny,nz,squeezeResult)
    % from  an array of 2D linear indices relative to a matrix of size [nx,ny], 
    % generates an array of linear indices of all the voxels of the 3D array
    % (size [nx ny nz]) which 2D projection is a member of the list of 2D
    % linear indices.
    % if the input has size [ni nj], the output has size [ni nj nz]
    % if the option squeezeResult has been selected, the dimensions of zie
    % 1 will be squeezed to return the lowest possible dimension array.

    % map idx2D array to a vector
    s = size(idx2D);
    idx2D  = idx2D(1:numel(idx2D));

    % Convert to subscript indices (i,j)
    [i, j] = ind2sub([nx, ny], idx2D);
    
    % Expand to all k values
    k = reshape(1:nz, 1, nz);  % [1 x nz]
    i = reshape(i, [], 1);     % [N x 1]
    j = reshape(j, [], 1);     % [N x 1]
    
    % Broadcast and get all (i,j,k)
    I = repmat(i, 1, nz);  % [N x nz]
    J = repmat(j, 1, nz);  % [N x nz]
    K = repmat(k, length(i), 1);  % [N x nz]
    
    % Convert back to linear indices into [nx ny nz]
    idx3D = sub2ind([nx ny nz], I(:), J(:), K(:));
    
    % recover original array format of the input
    idx3D = reshape(idx3D,[s,nz]);
    
    % remove useless dimensions if option has been selected
    if squeezeResult
        idx3D = squeeze(idx3D);
    end
end

%%
function paddedImg = padWithObjectValues(img, mask, maskVal)
    
    % check that size and dimensions of img and mask are compatible 
    [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask);
    if errorFlag == 1
        paddedImg = [];
        return
    end

    % Ensure mask is logical
    img = double(img);
    mask = (mask == maskVal);

    % Compute distance transform and index of nearest object pixel
    [~, idx] = bwdist(mask);
    
    % if mask is 2D and image is 3D, convert the 2D closest pixel in the
    % mask into its 3D counterpart
    if ndims(img) == 3 && ismatrix(mask)
        squeezeResult = 0;
        [nx,ny,nz] = size(img);
        idx = expand2DIndicesTo3D(idx,nx,ny,nz,squeezeResult);
    end

    % Initialize result
    paddedImg = img;

    % Fill background pixels with values from nearest object pixels
    paddedImg(~mask) = img(idx(~mask));
end

%%
function [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask)
    % image and mask must be 2D or 3D
    % nx ny must be the same.
    % if img is 2D, mask needs to be 2D
    % if img is 3D, mask can either be 2D or have the same nz as img.
    ndImg = ndims(img);
    ndMask = ndims(mask);
    errorFlag = 0;

    if ~ismember(ndImg,[2,3])
        disp(['Error: Image is ',num2str(ndImg),'D; needs to be 2D or 3D.']);
        errorFlag = 1;
        return
    end
    if ~ismember(ndMask,[2,3])
        disp(['Error: Mask is ',num2str(ndMask),'D; needs to be 2D or 3D.']);
        errorFlag = 1;
        return
    end

    if size(img,1) ~= size(mask,1) || size(img,2) ~= size(mask,2)
        disp(['Error: Image 2D size ',num2str([size(img,1),size(img,2)]),' and mask 2D size ',...
            num2str([size(mask,1),size(mask,2)]),' should be equal;',...
            'cannot find mask indices in img.']);
        errorFlag = 1;
        return
    end

    if ndImg == 2 && ndMask == 3
        disp(['Error: Image dimensions ',num2str(ndImg),' and mask dimensions ',...
            num2str(ndMask),' are incompatible;',...
            'cannot find mask indices in img.']);
        errorFlag = 1;
        return
    end

    if ndImg == 3 && ndMask == 3
        if size(img,3) ~= size(mask,3) 
            if size(mask,3) == 1
                % convert mask to a bona fide 2D matrix if the third
                % dimension has size 1.
                mask = squeeze(mask);
            else
                disp(['Error: Image z size ',num2str(size(img,3)),' and mask z size ',...
                    num2str(size(mask,3)),' should be equal;',...
                    'cannot find mask indices in img.']);
                errorFlag = 1;
                return
            end
        end
    end

end

%%
function smooth = smoothImg(img, numDim, filterHi, filterLo, psfSigma_xy)
    if numDim == 3 || numDim == 2    
        
        %factors 1.5 is a heuristic attempt to reproduce results
        %from the Fourier filter given the same parameters
        smooth = 1.5*DOGfilter(img, filterHi, filterLo * psfSigma_xy, [],[]);
    else
        disp('data type error: smoothing only possible on 2d images or 3d stacks');
        smooth = 0;
        return;
    end
end

%%
function [m,s,m2,s2] = getMeanStdInMask(img,mask,maskVal)
    % uses an ROI of ID value maskVal in the reference mask array mask to
    % compute summary statistics in the matching pixels/voxels of img.
    % Outputs the mean (m), std (s), 
    % the robust estimate of the mode (m2), 
    % and the pseudo std (s2) which is obtained as m2 - median( pixels with values below m2)
    % the function handles img and mask having the same size, 
    % or img being 3D and mask being 2D, in which case all img voxels that project
    % vertically onto the ROI in the mask are included in the summary
    % statistics.

    % check dimensions compatibility
    [errorFlag,mask] = checkImgMaskSizeDimMismatch(img,mask);
    if errorFlag ==1
        m = NaN; s = NaN;m2 = NaN; s2 = NaN;
        return
    end
    
    % number of bins used to compute the mode 
    % (these bins will map the 20th 80th percentile interval)
    nBins = 20;
    
    idx = getIdxInImgMatchingMaskValue(img,mask,maskVal);
    img = img(idx);
    m = mean(img(:));
    s = std(img(:));
    m2 = estimateModeRobust(img(:),nBins);
    s2 = m2 - median(img(img <= m2));

end

%%
function xMode = estimateModeRobust(d,nBins)
    % estimates the mode of a distribution of data d
    % by generating a histogram and fitting the max of the histogram to a
    % quadratic function to foind the max

    % INPUT
    % d is the data
    % nBins is the number of bins within the 20th-80th percentile interval
    % OUTPUT
    % xMode is the estimated value of the mode

    % generates a binning grid so that there are nBins between the 20th and
    % 80th percentile of the dataset data.
    binEdges = computeRobustHistBins(d,nBins); 
    [n,x] = hist(d,binEdges);
    
    % estimate the distribution mode
    [~,idxMode] = max(n); % this extracts the index of the maximum bin
    winSize = round(nBins/4);
    idxw = max(1,idxMode-winSize) : min(numel(x),idxMode+winSize); 

    % use a quadratic approximation of the histogram datapoints around the max
    % to estimate the distribution mode xModeFit
    nw = n(idxw);
    xw = x(idxw);
    [xMode,~] = findQuadMax2(nw',xw');
end

%% find approximate x,y coordinate of local max
function [xMax, yMax] = findQuadMax2(y,x)
    % given a series of x,y values assumed to bracket a local max
    % fits y(x) to a second degree polynomial y = ax^2 + bx +c
    % then outputs the extremum value xMax, i.e. x such that dy/dx = 0, i.e. x = -b/(2a)
    
    % Ensure column vectors
    x = x(:);
    y = y(:);

    % Precompute required sums
    S0 = length(x);
    S1 = sum(x);
    S2 = sum(x.^2);
    S3 = sum(x.^3);
    S4 = sum(x.^4);
    T0 = sum(y);
    T1 = sum(x .* y);
    T2 = sum((x.^2) .* y);

    % Coefficient matrix
    A = [S4, S3, S2;
         S3, S2, S1;
         S2, S1, S0];

    % Right-hand side vector
    B = [T2; T1; T0];

    % Solve the system
    warning('off', 'MATLAB:nearlySingularMatrix');
    coeffs = A \ B;
    warning('on', 'MATLAB:nearlySingularMatrix');
    % Extract coefficients
    a = coeffs(1);
    b = coeffs(2);
    c = coeffs(3);
    
    xMax = -b/(2*a);
    yMax = a*xMax^2 + b*xMax + c;
end

%%
function binEdges = computeRobustHistBins(data,nBins)
    % generates a binning grid so that there are nBins between the 20th and
    % 80th percentile of the dataset data.

    % outputs the edges of the bins, which encompass the entire dataset
    % (this means the number of bins in the output will be greater than nBins).

    % This approach ensure robust sampling of the mode of the distribution 
    % even for long-tailed distributions.

    data = data(~isinf(data));

    % Filter data within the 20-80 percentile range
    p20 = prctile(data, 20) 
    p80 = prctile(data, 80) 
    
    % split that range into nBins
    binSize = (p80-p20)/nBins
    
    % Define bin edges (n equal-sized bins)
    binEdges = min(data):binSize:max(data);
end

%%
function mList = collectMaskIds(mask,eliminateBackgroundROIs,backgroundID)
    mList = unique(mask(:));
    
    % eliminate background spots if needed
    if eliminateBackgroundROIs
        if isnan(backgroundID)
            mList = setdiff(mList,min(mList(:)));
        else
            mList = setdiff(mList,backgroundID);
        end
    end
end

%%
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
