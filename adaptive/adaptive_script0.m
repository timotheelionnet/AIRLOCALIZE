% adaptive script
% paths
imgName = 'examples/adaptive_examples/MFGTMP_250324110001_B02f00d1.TIFF';
maskName = 'examples/adaptive_examples/MFGTMP_250324110001_B02f00d1_cp_masks.tif';

% options

% if selected, eliminates spots in the background region
eliminateBackgroundSpots = 1;

% pixel value that marks the backgorund spots. If set to NaN, the region
% with minimal pixel value will be eliminated.
backgroundID = 0;

% factor that will be multiplied to the std of the mask
threshFactor = 2;

% adding 10 pixels around each mask in the cropped image
paddingSize = 0;

%% load image and mask
[img,nImgSlices] = timtiffread(imgName);
[mask,nMaskSlices] = timtiffread(maskName);

% collect masks IDs
mList = collectMaskIds(mask,eliminateBackgroundSpots,backgroundID);
nm = numel(mList);

%% initialize Airlocalize structures
alData = airLocalizeData();
alData.curFile = imgName;
alData.img = img; % holds the raw image/stack data of the current file
alData.isMovie = 0; 

params = airLocalizeParams();
params.reset;
params.threshUnits = 'absolute'; % units of the threshold - absolute or SD
params.minDistBetweenSpots = 2; %: minimum allowed distance between local maxima
params.numDim = ndims(img); %: the number of dimensions of the image (2 or 3)
params.outputSpotsImage = 1;
params.psfSigma = 1.3;
params.fitMethod = '2DGaussianMask';
verbose = 1;
%% loop through masks
maskStats = zeros(nm,4);
loc = [];
for i=1:nm
%for i=190:190
    % get mean & std of intensity over the mask region
    [m,s,m2,s2] = getMeanStdInMask(img,mask,mList(i));
    maskStats(i,:) = [m,s,m2,s2];

    % crop mask + padding
    [croppedImg, dx, dy, dz] = cropMask(img, mask, mList(i), paddingSize);
    
    % compute threshold based on mask intensity
    %params.threshLevel = m2 + threshFactor*s2;
    params.threshLevel = threshFactor*s2;

    % run gaussian localization on cropped mask
    alData.img = croppedImg; % current img
    alData.smooth = smooth_image_and_subtract_background7(alData.img,params);
    %alData.smooth = alData.img;

    spotCandidates = find_isolated_maxima_clean3(...
            alData,params,verbose);   

    % detection/quantification
    [tmpLoc,locVars] = run_gaussian_fit_on_all_spots_in_image3(...
        spotCandidates,alData,params);
    tmpLoc(:,ndims(img)+3) = mList(i); % updating the "frame" column to the mask ID
    
    % map spot coordinates in cropped image back to entire image.
    cropOffset = [dx, dy, dz] - 1; % offset starts at zero
    tmpLoc(:,1:ndims(img)) = tmpLoc(:,1:ndims(img)) ...
        + repmat(cropOffset,size(tmpLoc,1),1);

    % eliminate spots that do not fall inside the mask
    if ndims(img) == 2
        idx = sub2ind( size(img), ceil(tmpLoc(:,1)), ceil(tmpLoc(:,2)) );
    else
        idx = sub2ind( size(img), ceil(tmpLoc(:,1)), ceil(tmpLoc(:,2)), ...
            ceil(tmpLoc(:,3)) );
    end
    isInMask = mask(idx) == mList(i);
    nDelete = size(tmpLoc,1) - sum(isInMask);
    disp(['Deleted ',num2str(nDelete),' spots out of mask']);
    tmpLoc = tmpLoc(isInMask,:);
    
    % compile spots into full list 
    loc = [loc;tmpLoc]; 

end

%%
% save spot coordinates and detection parameters to text file
params.saveLocAndPar(loc,locVars,alData);

% generate tiff image
if params.outputSpotsImage
    disp('  saving spots image...');
    alData.img = img;
    params.saveSpotImg(loc,locVars,alData);
end

%% functions

%%
function [m,s,m2,s2] = getMeanStdInMask(img,mask,maskVal)
    % outputs the mean (m), std (s), 
    % the robust estimate of the mode (m2), 
    % and the pseudo std (s2) which is obtained as m2 - median( pixels with values below m2)
    
    nBins = 20;
    
    img = img(mask == maskVal);
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
    p20 = prctile(data, 20); 
    p80 = prctile(data, 80); 
    
    % split that range into nBins
    binSize = (p80-p20)/nBins;
    
    % Define bin edges (n equal-sized bins)
    binEdges = min(data):binSize:max(data);
end

%%
function mList = collectMaskIds(mask,eliminateBackgroundSpots,backgroundID)
    mList = unique(mask(:));
    
    % eliminate background spots if needed
    if eliminateBackgroundSpots
        if isnan(backgroundID)
            mList = setdiff(mList,min(mList));
        else
            mList = setdiff(mList,backgroundID);
        end
    end
end

%%
function [croppedImg, dx, dy, dz] = cropMask(img, mask, mID, paddingSize)
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

    is3D = ndims(img) == 3;
    if size(img,1) ~= size(mask,1) || size(img,2) ~= size(mask,2)
        error('Mask must be the same size as image');
    end
    if is3D
        if ndims(mask) == 2
            % If mask is 2D but img is 3D, replicate along z
            mask = repmat(mask, [1, 1, size(img, 3)]);
        elseif size(img,3) ~= size(mask,3) 
            error('Mask must be 2D or the same size as img');
        end
    end

    % Logical mask of the selected region
    region = (mask == mID);

    % Find bounding box
    if is3D
        [x, y, z] = ind2sub(size(region), find(region));
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
