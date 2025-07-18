
function [loc,locVars,locParams] = localizeImgPairAdaptive(img,mask,varargin)
    % from a pair of a data image and a matching mask image, localizes
    % spots across all masks of the image, using an adaptive threshold that
    % is scaled to the intensity std within each mask.
    
    % the value of the threshold is set through name-value arguments, either 
        % 'params' (airlocalizeParams object which contains the threshLevel property)
        % or directly via 'threshLevel' argument (recommended value 6 in 2D, 12 is
        % 3D).
    % the size of the psf is set through name-value arguments, either 
        % 'params' (airlocalizeParams object which contains the psfSigma property)
        % or directly via 'psfSigma' argument.

    % INPUT
    % img: data image (can be a z-stack)
    % mask: mask image (can be a z-stack or a 2D image if
        % the data image is a 3D stack)
    
    % optional arguments (name-value pair formatted)
    % 'params': airLocalizeParams object that fitting options will be read
        % from.
    % 'psfSigma': in 2D, sigma_xy; in 3D [sigma_xy, sigma_z] in pixel units;
        % the width of the gaussian fitting the point spread function in 2D
        % or 3D depending on the dimensionality of the image. 
        % !!! THIS PARAMETER IS REQUIRED IF PARAMS IS NOT POPULATED !!!
        % if both params and psfSigma arguments are input and psfSigma isnt
        % empty, psfSigma overrides the value of the psfSigma property in
        % params.
    % 'threshLevel': the factor used to determine the intensity threshold of
        % the regions that are spot candidates. the threshold is set as
        % thresh = max(0,mode) + threshFactor * s; where mode is the
        % estimate of the mode of the intensity in the smoothed image of
        % each object; and s an estimate of the std of the intensisty in
        % the smoothed image of each object.
        % !!! THIS PARAMETER IS REQUIRED IF PARAMS IS NOT POPULATED !!!
        % if both params and threshevels arguments are input and threshLevel isnt
        % empty, threshLevel overrides the value of the threshLevel property in
        % params.
    % 'eliminateBackgroundSpots': set to 0 to eliminate background objects
        % (default is 1)
    % 'backgroundID': ID of the background objects
        % (default is 0)
    % 'paddingSize': size in voxels of the padding around each object used 
        % before cropping. Default (recommended): 0    
    % 'verbose': sets the verbose setting (0 or 1) of the airlocalize 
        % detection function (default 1).
    
    % OUTPUT
    % loc: npts * (ndims + 3) array listing the position and intensities of the spots, with
        % the format x, y, (z if 3D), Intensity, Intensity_Residuals, Object_ID
    % locVars: cell array holding the names of the variables in the loc array 
    % locParams
    
    % collect optional arguments
    p = inputParser();
    addParameter(p,'locParams',[]); 
    addParameter(p,'psfSigma',[]);
    addParameter(p,'threshLevel',[]);  
    addParameter(p,'eliminateBackgroundSpots',1);    
    addParameter(p,'backgroundID',0);    
    addParameter(p,'paddingSize',0);    
    addParameter(p,'verbose',1);    
    
    parse(p,varargin{:});
    paddingSize = p.Results.paddingSize;
    locParams = p.Results.locParams;
    psfSigma = p.Results.psfSigma;
    threshLevel = p.Results.threshLevel;
    eliminateBackgroundSpots = p.Results.eliminateBackgroundSpots;
    backgroundID = p.Results.backgroundID;
    verbose = p.Results.verbose;

    if isempty(locParams) && isempty(psfSigma)
        error('You need to fill either locParams or psfSigma');
    end
    if isempty(locParams) && isempty(threshLevel)
        error('You need to fill either locParams or threshLevel');
    end

    % initialize the Airlocalize structure that will hold the localization
    % settings
    nd = ndims(img);
    if isempty(locParams)  
        locParams = airLocalizeParams();
        locParams.reset;
        locParams.numDim = nd; %: the number of dimensions of the image (2 or 3)
        locParams.threshUnits = 'absolute'; % units of the threshold - absolute or SD
        if nd == 2
            locParams.fitMethod = '2DGaussianMask';
        else
            locParams.fitMethod = '3DMaskFull';
        end
    else
        % create 'deep copy' of the object 
        % (since locParams is a handle and we don't want to modify the original)
        locParams = copy(locParams);
    end
    if ~isempty(psfSigma)
        locParams.psfSigma = psfSigma;
    end
    if ~isempty(threshLevel)
        locParams.threshLevel = threshLevel;
    end
    % safeguarding the threshLevel value now because we will update 
    % the threshLevel property of locParams for each object during the loop
    threshLevel = locParams.threshLevel;

    % initialize Airlocalize structure that will hold the image data
    alData = airLocalizeData();
    alData.img = img; % holds the raw image/stack data of the current file
    alData.isMovie = 0; 

    % replicate the mask in the third dim if it is two-dimensional and the
    % image is a stack.
    if nd == 3 && ismatrix(mask)
        mask = repmat(mask,1,1,size(img,3));
    end
    
    % collect list of masks IDs
    mList = collectMaskIds(mask,eliminateBackgroundSpots,backgroundID);
    nm = numel(mList);
    
    %% loop through masks
    maskStats = zeros(nm,4);
    loc = [];
    for i=1:nm
        % crop mask + padding
        [croppedImg, dx, dy, dz] = cropMask(img, mask, mList(i), paddingSize);
        croppedMask = cropMask(mask, mask, mList(i), paddingSize);
        
        % run gaussian localization on cropped mask
        alData.img = croppedImg; % current img
        alData.smooth = ...
            smooth_image_and_subtract_background7(alData.img,locParams);
        
        % get mean & std of intensity over the mask region
        [m,s,m2,s2] = getMeanStdInMask(alData.smooth,croppedMask,mList(i));
        maskStats(i,:) = [m,s,m2,s2];
        
        % compute threshold based on mask intensity - using the average of
        % the two "Standard deviations" bc it is a bit more robust
        locParams.threshLevel = max(0,m2) + threshLevel*(s2+s)/2;
    
        spotCandidates = find_isolated_maxima_clean3(...
                alData,locParams,verbose);   
    
        % detection/quantification
        [tmpLoc,locVars] = run_gaussian_fit_on_all_spots_in_image3(...
            spotCandidates,alData,locParams);
        tmpLoc(:,nd+3) = mList(i); % updating the "frame" column to the mask ID
        
        % map spot coordinates in cropped image back to entire image.
        cropOffset = [dx, dy, dz] - 1; % offset starts at zero
        tmpLoc(:,1:nd) = tmpLoc(:,1:nd) ...
            + repmat(cropOffset,size(tmpLoc,1),1);
    
        % eliminate spots that do not fall inside the mask
        if nd == 2
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

    % update ROI_ID as the column of the loc table normally used for frame numbers
    idx = ismember(locVars,'image_number');
    if sum(idx) ~= 0
        locVars{idx} = 'ROI_ID';
    end
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

    nd = ndims(img);

    is3D = nd == 3;
    if size(img,1) ~= size(mask,1) || size(img,2) ~= size(mask,2)
        error('Mask must be the same size as image');
    end
    if is3D
        if ismatrix(mask)
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
