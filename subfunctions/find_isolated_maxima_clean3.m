function maxima = find_isolated_maxima_clean3(alData,params,verbose)
% finds local maxima above the threshold in the smoothed image

% alData is an airlocalize object that needs to contain 
    % the smoothed image alData.smooth
    % the original image alData.img (used if threshold is in units of SD)
    % the field isMovie that instructs whether to treat just the current
    % frame
    % the curFrame field that instruct which frame to process.
% params is an airlocalizeParams object that needs to contain the following properties: 
    % threshLevel: value of the threshold
    % threshUnits: units of the threshold - absolute or SD
    % minDistBetweenSpots: minimum allowed distance between local maxima
    % numDim: the number of dimensions of the image (2 or 3)

% outputs maxima: a list of local maxima pixel coordinates
% [ x y I ] in 2D
% [ x y z I ] in 3D

%% retrieve threshold level 

if strcmp(params.threshUnits,'absolute')
    threshInt = params.threshLevel;
    if verbose
        disp(['  threshold value is ' num2str(threshInt) ' in absolute units']);
    end
elseif strcmp(params.threshUnits,'SD')
    if ~alData.isMovie
        threshInt = params.threshLevel * std(alData.img,0,'all');
    else
        threshInt = params.threshLevel ...
            *std(alData.img(:,:,alData.curFrame),0,'all');
    end
    if verbose
        disp(['  threshold value is ', num2str(params.threshLevel),...
        ' in SD units, i.e. ', num2str(threshInt),...
        ' in absolute units for current frame']); 
    end
elseif strcmp(params.threshUnits,'legacySD')
    if ~alData.isMovie
        threshInt = params.threshLevel * std(alData.smooth,0,'all');
    else
        threshInt = params.threshLevel ...
            *std(alData.smooth(:,:,alData.curFrame),0,'all');
    end
    if verbose
        disp(['  threshold value is ', num2str(params.threshLevel),...
        ' in (legacy) SD units, i.e. ', num2str(threshInt),...
        ' in absolute units for current frame']); 
    end
end

%% find local maxima
minDist = round(params.minDistBetweenSpots );
if minDist >0
    %find local maxima within ROI pixels in each dimension
    if ndims(alData.smooth) == 3 && ~alData.isMovie
        %generate local window for dilation operation
        dilwin = ones(2*minDist+1,2*minDist+1,2*minDist+1);
        dilwin(minDist+1,minDist+1,minDist+1) = 0;

    else
        %generate local window for dilation operation
        dilwin = ones(2*minDist+1,2*minDist+1);
        dilwin(minDist+1,minDist+1) = 0;
    end
    % imdilate result is an image in which each pixel takes the value
    % of the highest pixel (except itself) in its neighborhood
    if ~alData.isMovie
        localmax = alData.smooth >= imdilate(alData.smooth,dilwin);
    else
        localmax = alData.smooth(:,:,alData.curFrame) ...
            >= imdilate(alData.smooth,dilwin);
    end
else
    if ~alData.isMovie
        localmax = ones(size(alData.smooth));
    else
        localmax = ones(size(alData.smooth(:,:,alData.curFrame)));
    end
end

% enforce that local maxima intensity is above threshold 
localmax = localmax.*(alData.smooth > threshInt);

%% store maxima as a list of coordinates / intensity
x = find(localmax);
clear('localmax');
Int = alData.img(x);
if ndims(alData.img) == 3
    [x,y,z] = ind2sub(size(alData.img),x);
    maxima = [x,y,z,Int];
else
    [x,y] = ind2sub(size(alData.img),x);
    maxima = [x,y,Int];
end

if isempty(maxima)
    if verbose
        disp('  predetected no spots;'); 
    end
    return
end

% ordering the maxima by descending intensity value
maxima = sortrows(maxima,params.numDim+1,'descend');

%truncating the array if it has more points than allowed
if size(maxima,1) > params.maxSpots
    maxima = maxima(1:params.maxSpots,:);
    if verbose
        disp( ['  predetected ',num2str(size(maxima,1)),...
            ' spots, reducing their number to the max allowed number = ',...
            num2str(params.maxSpots)]);
    end
else
    if verbose
        disp(['  predetected ',num2str(size(maxima,1)),' spots;']);
    end
end
    
end