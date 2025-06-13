
function loc = localizeFilePairAdaptive(imgName,maskName,psfSigma,threshLevel,...
                    varargin)
    % from a pair of a data image and a matching mask image, localizes
    % spots across all masks of the image, using an adaptive threshold that
    % is scaled to the intensity std within each mask.

    % INPUT
    % imgName: file path of the data image (can be a z-stack)
    % maskName: file path of the mask image (can be a z-stack or a 2D image if
        % the data image is a 3D stack)
    % psfSigma: in 2D, sigma_xy; in 3D [sigma_xy, sigma_z] in pixel units;
        % the width of the gaussian fitting the point spread function in 2D
        % or 3D depending on the dimensionality of the image.
    % threshLevel: the factor used to determine the intensity threshold of
        % the regions that are spot candidates. the threshold is set as
        % thresh = max(0,mode) + threshFactor * s; where mode is the
        % estimate of the mode of the intensity in the smoothed image of
        % each object; and s an estimate of the std of the intensisty in
        % the smoothed image of each object.
    
    % optional (name-value pair formatted)
    % 'saveDirName': path to the saving directory (default is the folder
        % holding the image)
    % 'eliminateBackgroundSpots': set to 0 to eliminate background objects
        % (default is 1)
    % 'backgroundID': ID of the background objects
        % (default is 0)
    % 'generateSpotsImage': set to 1 to generate an image of the spots,0
        % otherwise.

    % OUTPUT
    % loc: table listing the position and intensities of the spots, with
    % the format x,y,(z),Intensity,Intensity_Residuals,Object_ID
    
    % collect optional arguments
    p = inputParser();
    addParameter(p,'saveDirName','');    
    addParameter(p,'eliminateBackgroundSpots',1);    
    addParameter(p,'backgroundID',0);    
    addParameter(p,'generateSpotsImage',1);

    parse(p,varargin{:});
    saveDirName = p.Results.saveDirName;
    eliminateBackgroundSpots = p.Results.eliminateBackgroundSpots;
    backgroundID = p.Results.backgroundID;
    generateSpotsImage = p.Results.generateSpotsImage;
    
    % load image and mask from files
    img = timtiffread(imgName);
    mask = timtiffread(maskName);
    
    [loc,locVars,params] = localizeImgPairAdaptive(img,mask,...
                    'eliminateBackgroundSpots',eliminateBackgroundSpots,...
                    'backgroundID',backgroundID,...
                    'psfSigma',psfSigma,'threshLevel',threshLevel);
    
    % save spot coordinates and detection parameters to text file
    if ~isempty(saveDirName)
        params.saveDirName = saveDirName;
    else
        params.saveDirName = fileparts(imgName);
    end
    alData = airLocalizeData();
    alData.img = img;
    alData.curFile = imgName;
    params.saveLocAndPar(loc,locVars,alData);
    
    % generate tiff image
    params.outputSpotsImage = generateSpotsImage;
    if params.outputSpotsImage
        disp('  saving spots image...');
        alData.img = img;
        params.saveSpotImg(loc,locVars,alData);
    end
end

