function loc = AIRLOCALIZE(varargin)
% GUI-less operation: 
    % AIRLOCALIZE(parametersFilename) 
    % loads a text file holding all parameters 
    % needed for analysis

% GUI mode: lets the user interactively set all parameters
    % AIRLOCALIZE()  (a pop up will appear for file loading)
    % GUI mode is compatible with batch analysis.
    
% AIRLOCALIZE_version = 1.6;  

addpath('subfunctions');
addpath('subfunctions/iniconfig');

persistent params;
loc = [];
%% collect parameters
if nargin == 0  
    % if no input selected, enter parameters via GUI
    [params,alData] = set_detection_parameters3(params); 
    if strcmp(params.fileProcessingMode,'cancel'), return;  end 
    
elseif nargin == 1
    % read parameters from config file
    disp('reading config file...');
    parametersFilename = varargin{1};
    if ~isa(params,'airlocalizeParams') 
        params = airLocalizeParams();
        params.reset();
    end
    [params,status] = params.read_params_from_file2(parametersFilename);
    if strcmp(status,'cancel'), return;  end 
    
    % get list of files to analyze
    disp('collecting list of files to analyze...');
    alData = airLocalizeData();
    alData.setFListFromParams(params);
    if isempty(alData.getFList)
        disp('could not find files to analyze.');
        return
    else
        disp(['Found ',num2str(numel(alData.getFList)),' files to analyze.']);
    end
    
    % if empty, default saving dir to the dir holding the first image in the
    % list
    if isempty(params.saveDirName)
        params.saveDirName = fileparts(alData.fList{1});
    end

    % set movie mode
    if contains(params.fileProcessingMode,'movie') ...
            || contains(params.fileProcessingMode,'Movie')
        alData.isMovie = 1;
    end
    
    % figure out data dimensionality on first file of the list.
    disp('setting data dimensionality...');
    params.getNumDim(alData.getFList);
    if strcmp(params.fileProcessingMode,'cancel'), return;  end
    
    % check parameters consistency
    params.checkParamsConsistency(alData,'numDim');
    
    % set PSF width manually (optional)
    alData.setFileIdx(1);
    if params.setPsfInteractively == 1
        disp('setting PSF size interactively...');
        [params,alData] = set_psf_size_manually(params,alData);
    end
    if strcmp(params.fileProcessingMode,'cancel'), return; end

    % set Int threshold manually (optional) 
    alData.setFileIdx(1);
    if params.setThreshInteractively  == 1 
       disp('setting threshold interactively...');
       params = set_threshold_manually(params,alData); 
    end
    if strcmp(params.fileProcessingMode,'cancel'), return; end
end

%% analyze file(s)

% create saving dir if necessary
if ~isempty(params.saveDirName)
    if ~exist(params.saveDirName,'dir')
        disp(['Creating saving directory ',params.saveDirName,' ...']);
        mkdir(params.saveDirName);
    end
end

for i=1:numel(alData.fList)
    % update the current file in the list
    alData.setFileIdx(i);
    disp(['analyzing file: ',alData.curFile,'...']);
    
    % retrieve current image and smooth it
    if i == 1
        overwrite = 0;
    else
        overwrite = 1;
    end
    alData.retrieveImg(overwrite);
    
    % verify that dimensionality agrees with settings
    skipFile = 0;
    nd  = ndims(alData.img);
    if nd == 2
        if params.numDim == 3
            disp(['  Current image has ',num2str(nd),' dimensions but expecting 3D data; skipping file']);
            skipFile = 1;
        end
        if alData.isMovie
            disp('  Warning: file expected to be a movie but 2D image was found.');
        end
    elseif nd == 3
        if params.numDim == 2 && ~alData.isMovie
            disp(['  Current image has ',num2str(nd),' dimensions but expecting 2D data; skipping file']);
            skipFile = 1;
        end
    else
        skipFile = 1;
        disp(['  Current image has ',num2str(nd),' dimensions; cannot analyze; skipping file']);
    end
    
    if skipFile
        continue;
    end
    
    alData.retrieveSmoothedImg(params,overwrite);
    % go through all frames (if file is an image, just one "frame")
    loc = [];
    locVars = [];
    for j = 1 : alData.nFrames
        if alData.isMovie
            disp(['  analyzing frame ',num2str(j),'...']);
        end
        if alData.nFrames > 100
            if mod(j,100) == 0
                disp(['  analyzing frame ',num2str(j),'...']);
            end
        end
        
        % update current frame
        alData.setCurFrame(j);
        
        % predetection of local maxima
        verbose = 1;
        spotCandidates = find_isolated_maxima_clean3(...
            alData,params,verbose);   

        % detection/quantification
        [tmpLoc,locVars] = run_gaussian_fit_on_all_spots_in_image3(...
        spotCandidates,alData,params);
        
        loc = [loc; tmpLoc];
    end
    
    % save spot coordinates and detection parameters to text file
    params.saveLocAndPar(loc,locVars,alData);
    
    if params.outputSpotsImage
        disp('  saving spots image...');
        params.saveSpotImg(loc,locVars,alData);
    end
    
    disp('  done analyzing file.');
end
disp('done.');

end



