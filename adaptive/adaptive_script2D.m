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
imgName = 'examples/adaptive_examples/MFGTMP_250324110001_B02f01d1.TIFF';
maskName = 'examples/adaptive_examples/MFGTMP_250324110001_B02f01d1_cp_masks.tif';

%% settings
% key parameter: size of the psf in pixel units 
% should be scalar for 2D image data / should be 2D for 3D z-stacks ([sigma_xy, sigma_z])
psfSigma = 1.3;

% key parameter: factor that will be multiplied to the std of the mask
% (recommended: 6)
threshFactor = 6;

% options
% if = 1, eliminates spots in the background region; set to 0 to keep (default: 1)
eliminateBackgroundSpots = 1;

% pixel value that marks the backgorund spots. If set to NaN, the region
% with minimal pixel value will be eliminated (default: 0)
backgroundID = 0;

% size of the pixel padding around each mask in the cropped image (default:
% 0)
paddingSize = 0;

%% run the localization on the file pair
% loc: table listing the position and intensities of the spots, with
    % the format x,y,(z),Intensity,Intensity_Residuals,Object_ID
loc = localizeFilePairAdaptive(imgName,maskName,psfSigma,threshFactor,...
                    eliminateBackgroundSpots,backgroundID,paddingSize,...
                    'saveDirName','/Users/lionnt01/Documents/junk/');

