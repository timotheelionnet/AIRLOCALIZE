classdef airLocalizeParams < handle
    % object that manages access to Airlocalize parameters
    
    properties (GetAccess = 'public', SetAccess = 'public')
        %% [Files]  

        % path to data file name / folder name for batch analysis
        dataFileName;

        % path to saving directory
        saveDirName;

        % whether to process an individual image or files in batch; whether files are movies or not.    
        fileProcessingMode;

        % comma-separated string containing patterns that need to be present in the
        % full file path for it to be analyzed
        inclusionString;

        % comma-separated string containing patterns that if present in the full 
        % file path, will exclude the file from analysis
        exclusionString;

        % set to 1 to explore subdirectories of the parent folder (batch mode
        % only)
        recursive;

        %% [detectionSettings]
        psfSigma;
        threshUnits;
        threshLevel;
        maxSpots;
        fitMethod;

        %% [userInput]
        setPsfInteractively; 
        setThreshInteractively; 
        movieFrameUsedToSetParameters; 

        %% [outputImg]
        outputSpotsImage;
        autoSpotSize;
        spotSize;
        autoSpotIntensity;
        spotIntensity;

        %% [advanced]
        minDistBetweenSpots;
        filterLo;
        filterHi;
        fittedRegionSize; 
        tol; 
        bgRegionThickness;
        psfType;
        maxIterations; 
        bgCorrectionMode;
        useLegacyFormats;

        %% internal use
        numDim; % dimensionality is set internally by checking the size of the first image (getNumDim)

    end
    
    methods
        
        % constructor, all properties are empty
        function obj = airLocalizeParams()
        end
        
        % sets params to default
        function obj = reset(obj)
            d = obj.getSettingsDefaultAndPresets;
            f = fieldnames(d);
            for i=1:numel(f)
                obj.(f{i}) = d.(f{i});
            end
        end
        
        % read parameters to
        function [obj,status] = read_params_from_file2(obj,filename)
            % reads an ini config file (filename should be the path to it)
            % and updates AIRLOCALIZE settings
            % then
            % checks that all necessary parameters are present
            % checks that format of the parameters is correct
            % outputs the parameters structure.
            status = '';

            % checking that the file exists
            if ~exist(filename,'file')
                disp(['could not find config file: ',filename]);
                status = 'cancel';
                return
            end

            % read parameters from ini config file
            disp(['reading settings from config file: ',filename,' ...']);
            ini = IniConfig();
            ini.ReadFile(filename);

            sections = ini.GetSections();
            for i=1:numel(sections)
                keys = ini.GetKeys(sections{i});
                for j= 1:numel(keys)
                    obj.(keys{j}) = ini.GetValues(sections{i}, keys{j});
                end
            end

            disp(['checking file input settings: ',filename,' ...']);
            [obj,status] = obj.checkParametersFormat();
        end
        
        % duplicate object
        function dup = duplicate(obj)
            dup = airLocalizeParams();
            p = properties('airLocalizeParams');
            for i=1:numel(p)
                dup.(p{i}) = obj.(p{i});
            end
        end

        % save loc and par files corresponding to current image in alData
        function saveLocAndPar(obj,loc,locVars,alData)
            % generate file names
            [locFileName,parFileName] = generateLocParFileNames(obj,alData);

            % save loc file
            t = array2table(loc,'VariableNames',locVars);
            [~,~,e] = fileparts(locFileName);
            if strcmp(e,'.loc4')
                % save as tab-delimited file with headers
                writetable(t,locFileName,'Delimiter','\t','FileType','text');
            else
                % legacy format (no headers)
                save(locFileName,'loc','-ascii');
            end

            % save parameter file
            [~,~,e] = fileparts(parFileName);
            if strcmp(e,'.par4')
                % replace dataFileName (which can be the input folder name 
                % if batch with current data file name
                p2 = obj.duplicate();
                p2.dataFileName = alData.curFile;
                
                % save as ini file
                overwrite = 1;
                p2.saveParAsIni(parFileName,overwrite);
            else
                % use legacy format
                params2 = obj.legacyParamsConversion([],alData,'toLegacy');
                save_structure_as_text2(params2,parFileName);
            end
        end
        
        function saveParAsIni(obj,parFileName,overwrite)
            if overwrite && exist(parFileName,'file')
                delete(parFileName);
            end
            % convert params object to ini
            ini = IniConfig();
            
            % collect the keys to save the objects under
            [~,~,~,~,section] =  obj.getSettingsDefaultAndPresets();
            f = fieldnames(section);
            for i=1:numel(f)
                curSection = section.(f{i});
                if ~strcmp(curSection,'internal')
                    if ~ini.IsSections(curSection)
                        ini.AddSections({curSection});
                    end
                    if ~ini.IsKeys(curSection, f{i})
                        ini.AddKeys(curSection,f(i));
                    end
                
                    ini.SetValues(curSection, f{i}, obj.(f{i}));
                end
            end
            
            % save ini as iniConfig file
            ini.WriteFile(parFileName);
        end
        
        function saveSpotImg(obj,loc,locVars,alData)
            % generate output file name
            [~,~,spotsFileName] = obj.generateLocParFileNames(alData);
            
            % generate spots image
            spotsImg = generateSpotsImage(obj,loc,locVars,alData);
            
            % save spots image
            save_as_tiff(spotsImg,spotsFileName);
        end
        
        % generates the file names for saving par and loc files
        function [locFileName,parFileName,imgFileName] = ...
                generateLocParFileNames(obj,alData)
            if ~isprop(obj,'useLegacyFormats')
                defaults = obj.getSettingsDefaultAndPresets();
                obj.useLegacyFormats = defaults.useLegacyFormats;
            end

            [d,f,~] = fileparts(alData.curFile);

            % use saving directory if specified, 
            % otherwise analysis files will be saved next to original images
            if ~isempty(obj.saveDirName)
                d = obj.saveDirName;
            end

            if obj.useLegacyFormats
                locFileName = fullfile(d,[f,'.loc3']);
                parFileName = fullfile(d,[f,'.par']);
            else

                locFileName = fullfile(d,[f,'.loc4']);
                parFileName = fullfile(d,[f,'.par4']);
            end
            imgFileName = fullfile(d,[f,'_spots.tif']);
        end

        % checks that all needed entries in the parameters structure are present.    
        function [obj,status] = checkParametersFormat(obj)
            % checks that all needed entries in the parameters structure are present.
            % if absent, replaces by the default setting
            % checks that each entry has the right format (e.g. string vs numeric)
            % checks that each entry fulfills expected range of values (e.g. positive for a numeric,
            % or one of the preset modes for a string).

            status = '';

            % pull the expected field names for AIRLOCALIZE settings and their defaults
            [defaults,type,presets,needed] =  obj.getSettingsDefaultAndPresets();

            missingSettings = 0;
            fields = fieldnames(defaults);
            for i=1:numel(fields)
                % check if all expected parameters are present in p; 
                % set missing parameters to defaults when possible
                if ~isprop(obj,fields{i})
                    [obj,m] = obj.checkIfSettingIsRequiredAndSetToDefault(...
                        fields{i},defaults,needed,'missing');
                    missingSettings = max(m,missingSettings);
                else
                    % check that the parameter format extracted from p 
                    % matches the expected data type.
                    if ~isa(obj.(fields{i}),type.(fields{i}){1})
                        [obj,m] = obj.checkIfSettingIsRequiredAndSetToDefault(...
                            fields{i},defaults,needed,'type');
                        missingSettings = max(m,missingSettings);
                    else
                        % check that the parameter value is in the expected range
                        obj = obj.checkThatEntryTypeIsCorrect(...
                            fields{i},type,defaults,presets);
                    end
                end
            end
            if missingSettings
                status = 'failed';
                return
            else
                status = 'success';
            end
        end
        
        function [obj,status] = checkParamsConsistency(obj,alData,source)
        % checks that parameters stored in params are internally consistent
        % and consistent with data stored in airlocalizeData object

        % input: 
        % params:  airlocaliza params structure
        % alData: airlocalizeData object (or left empty if source is
            % 'detection_mode')
        % source can be either
            % 'detection_mode': if checkParamsConsistency is called after
                % input of the file Processing mode.
            % 'numDim': if checkParamsConsistency is called after
                % figuring out the dimensionality.
            % 'detection_parameters_interface': if checkParamsConsistency is 
                % called after setting the parameters
            % 'analysis': if checkParamsConsistency is called 
                % before analyzing one image.   

        % output:        
        % status: either 
            % 'success' if no inconsistency was found,
            % 'fixed' if an inconsistency was found but fixed, 
            % 'failed' if there was an inconsistency that can't be fixed.
        % params is the updated params object.

        % check first that all params are present and have the expected format.
        [obj,status] = obj.checkParametersFormat();
        if strcmp(status,'failed')
            status = 'failed; incorrect parameter format';
            return
        end

        paramsToCheck = {...
            'inclusionString',...
            'exclusionString',...
            'numDim',...
            'psfSigma',...
            'fitMethod',...
            'spotSize',...
            'psfType',...
            };

        for i=1:numel(paramsToCheck)
            % checks whether current paramsToCheck is compatible with the rest of
            % the settings
            status = ...
                obj.checkCompatibilitySingleSetting(...
                paramsToCheck{i},alData);

            % if a problem was flagged, takes action based on what source was
            % calling the function: either update params to make it compatible,
            % or trigger an error. status is set to success or failed based on
            % whether a remedy was found
            if contains(status,'failed')
                [obj,status] = obj.repairIncompatibleSetting(...
                    paramsToCheck{i},source);
            end
        end

        end

        function [obj,status] = repairIncompatibleSetting(obj,...
            paramToCheck,source)
            % checks whether incompatible parameters can be repaired, depending on the
            % source
            % source can be either
                % 'detection_mode': if parent function checkParamsConsistency 
                    % is called after input of the file Processing mode.
                % 'numDim': if parent function checkParamsConsistency is called after
                    % figuring out the dimensionality.
                % 'detection_parameters_interface': if parent function 
                    % checkParamsConsistency is called after setting the parameters
                % 'analysis': if parent function checkParamsConsistency is called 
                    % before analyzing one image.   
            %outputs:
            % status: either 
                % 'success' if no inconsistency was found,
                % 'fixed' if an inconsistency was found but fixed, 
                % 'failed' if there was an inconsistency that can't be fixed.
            % params is the updated params object.

            % get parameters defaults
            d =  obj.getSettingsDefaultAndPresets();

            switch paramToCheck
                case 'inclusionString'
                    % remove unnecessary inclusionString field when fileProcessingMode
                    % is set to single image.
                    obj.inclusionString ='';
                    status = 'fixed';
                    return
                case 'exclusionString'
                    % remove unnecessary inclusionString field when fileProcessingMode
                    % is set to single image.
                    obj.exclusionString ='';
                    status = 'fixed';
                    return
                case 'numDim'
                    % inconsistencies between data size and numdim cannot be fixed
                    switch source
                        case 'numDim'
                            status = 'failed';
                            return
                        case 'analysis'
                            status = 'failed';
                            return
                    end
                case 'psfSigma'
                    switch source
                        case 'numDim'
                            % at this stage of Airlocalize, values are 
                            % placeholders so replacing with the defaults here.
                            disp('Initializing psfSigma to default.');
                            obj.psfSigma = d.psfSigma;
                            if obj.numDim == 2
                                obj.psfSigma = obj.psfSigma(1);
                            end
                            status = 'fixed';
                            return
                        case 'analysis'
                            status = 'failed';
                            return
                    end
                case 'fitMethod'
                    switch source
                        case 'detection_parameters_interface'
                            if ismember(obj.fitMethod,...
                                    {'3dMaskFull','2dMaskOnLocalMaxProj','3DGaussianFit'}) ...
                                    && obj.numDim == 2
                                disp('setting fitMethod to default: 2DGaussianMask');
                                obj.fitMethod = '2DGaussianMask';
                                status = 'fixed';
                                return
                            elseif ismember(obj.fitMethod,...
                                    {'2DGaussianMask','2DGaussianFit'}) ...
                                    && obj.numDim == 3
                                disp(['setting fitMethod to default: ',d.fitMethod]);
                                obj.fitMethod = d.fitMethod;
                                status = 'fixed';
                                return
                            end
                        case 'analysis'
                            status = 'failed';
                            return
                    end

                case 'spotSize'
                    switch source
                        case 'numDim'
                            if obj.numDim == 3 && numel(obj.spotSize) > 2
                                obj.spotSize = obj.spotSize(1:2);
                                status = 'fixed';
                                return
                            elseif obj.numDim == 3 && numel(obj.spotSize) < 2
                                obj.spotSize = d.spotSize;
                                status = 'fixed';
                                return
                            elseif obj.numDim == 2 && numel(obj.spotSize) > 1
                                obj.spotSize = obj.spotSize(1);
                                status = 'fixed';
                                return
                            elseif obj.numDim == 2 && numel(obj.spotSize) < 1
                                obj.spotSize = d.spotSize(1);
                                status = 'fixed';
                                return
                            else
                                status = 'failed';
                                return
                            end
                        case 'detection_parameters_interface'
                            if obj.numDim == 3 && numel(obj.spotSize) > 2
                                obj.spotSize = obj.spotSize(1:2);
                                status = 'fixed';
                                return
                            elseif obj.numDim == 3 && numel(obj.spotSize) < 2
                                obj.spotSize = d.spotSize;
                                status = 'fixed';
                                return
                            elseif obj.numDim == 2 && numel(obj.spotSize) > 1
                                obj.spotSize = obj.spotSize(1);
                                status = 'fixed';
                                return
                            elseif obj.numDim == 2 && numel(obj.spotSize) < 1
                                obj.spotSize = d.spotSize(1);
                                status = 'fixed';
                                return
                            else
                                status = 'failed';
                                return
                            end
                        case 'analysis'
                            status = 'failed';
                            return
                    end
                case 'psfType'
                    switch source
                        case 'numDim'
                            disp('setting psfType to 2D default: integrated gaussian');
                            obj.psfType = 'integrated gaussian';
                            status = 'fixed';
                            return
                        case 'analysis'
                            status = 'failed';
                            return
                    end
            end
            status = 'success';
        end

        function status = checkCompatibilitySingleSetting(obj,...
            paramToCheck,alData)

            % checks whether field paramToCheck of params object is compatible with
            % other entries (including the data loaded in alData)
            % returns a status flag that can be either:
                % success: everything is compatible
                % 'failed' or 'failed; <some more detailed explanation>'
                % most failures trigger an explanantion message displayed in the
                % command line.

            switch paramToCheck
                case 'inclusionString'
                    if contains(obj.fileProcessingMode,'single')
                        if ~isempty(obj.inclusionString)
                            status = 'failed';
                            return
                        end
                    end
                case 'exclusionString'
                    if contains(obj.fileProcessingMode,'single')
                        if ~isempty(obj.exclusionString)
                            status = 'failed';
                            return
                        end
                    end
                case 'numDim'
                    % checking that any data loaded into alData.img has the right
                    % dimensionality
                    if ~isempty(alData)
                        if ~isempty(alData.img)
                            if obj.numDim ~= ndims(alData.img)
                                disp(['error: params set to 3D but loaded image is ',...
                                    num2str(ndims(alData.img))','D']);
                                status = 'failed; numDim incompatible with data.';
                                return
                            end 
                        end
                    end
                    % checking that any data loaded into alData.smooth has the right
                    % dimensionality
                    if ~isempty(alData)
                        if ~isempty(alData.smooth)
                            if obj.numDim ~= ndims(alData.smooth)
                                disp(['error: params set to 3D but loaded smoothed image is ',...
                                    num2str(ndims(alData.smooth))','D']);
                                status = 'failed; numDim incompatible with smoothed image.';
                                return
                            end 
                        end
                    end
                case 'psfSigma'
                    if numel(obj.psfSigma) ~= obj.numDim-1 
                        disp(['warning: psfSigma setting has size ',...
                                num2str(numel(obj.psfSigma))...
                                ,'. Expecting ',num2str(obj.numDim-1)]);
                        status = 'failed';
                        return
                    end
                case 'fitMethod'
                    if ismember(obj.fitMethod,...
                            {'3dMaskFull','2dMaskOnLocalMaxProj','3DGaussianFit'}) ...
                            && obj.numDim == 2
                        disp(['warning, fitMethod is set to ',obj.fitMethod,...
                            ' which is not compatible with 2D']);
                        status = 'failed';
                        return
                    elseif ismember(obj.fitMethod,...
                            {'2DGaussianMask','2DGaussianFit'}) && obj.numDim == 3
                        disp(['warning, fitMethod is set to ',obj.fitMethod,...
                            ' which is not compatible with 3D']);
                        status = 'failed';
                        return
                    end

                case 'spotSize'
                    if obj.numDim ~= numel(obj.spotSize) +1
                        status = 'failed';
                        return
                    end
                case 'psfType'
                    if obj.numDim == 2
                        if strcmp(obj.psfType, 'integrated gaussian std z')
                            disp(['warning: psfType is set to integrated gaussian std z',...
                                ' which is not compatible with 2D data']);
                            status = 'failed';
                            return
                        end
                    end
            end
            status = 'success';    
        end

        % checks whether property propName is required and if allowed
        % replaces it to default value
        function [obj,missingRequiredSetting] = ...
            checkIfSettingIsRequiredAndSetToDefault(obj,propName,...
            defaults,needed,criterion)

            switch criterion
                case 'missing'
                    str = 'could not find required setting in config file';
                case 'type'
                    str = 'input does not match expected format for setting';
            end

            missingRequiredSetting = 0;
            switch needed.(propName)
                case 0
                    % optional setting, set to default
                    obj.(propName) = defaults.(propName);
                case 1
                    % required setting, display error message
                    disp(['error: ',str,': ',propName]);
                    missingRequiredSetting = 1;
                case 2 
                    % setting only required in some circumstances, display
                    % warning and set to static default
                    disp(['warning: ',str,': ',propName,'; using default value.']);
                    obj.(propName) = defaults.(propName);
                case 3
                    % setting needed but can be set to default, display
                    % warning and set to static default
                    disp(['warning: ',str,': ',propName,'; using default value.']);
                    obj.(propName) = defaults.(propName);
                case 4    
                    % setting needed but can be set to dynamic default,
                    % setting to static default as place holder
                    disp(['warning: ',str,': ',propName,'; using default value.']);
                    obj.(propName) = defaults.(propName);
            end

        end

        % checks that property propName has the expected format
        function obj = checkThatEntryTypeIsCorrect(...
                obj,propName,type,defaults,presets)
            % checking that extra constraints (e.g. integer or positive) are
            % fulfilled
            if numel(type.(propName)) > 1
                for i=2:numel(type.(propName))
                    switch type.(propName){i}
                        case '>0'
                            if obj.(propName) <= 0

                                disp(['Warning: setting ',propName,...
                                    ' input value was set to ',...
                                    num2str(obj.(propName)),...
                                    ' but expected to be positive. Using absolute value instead.']);
                                obj.(propName) = abs(obj.(propName));
                            end
                        case 'integer'
                            if floor(obj.(propName)) ~= obj.(propName)

                                pOld = obj.(propName);
                                obj.(propName) = floor(obj.(propName));
                                if obj.(propName) == 0 && ismember('>0',type.(propName))
                                    obj.(propName) = 1;
                                end
                                disp(['Warning: setting ',propName,...
                                    ' input value was set to ',num2str(pOld),...
                                    ' but expected to be integer. Using ',...
                                    num2str(obj.(propName)),' instead.']);
                            end
                    end
                end 
            end

            % checking that if setting entry is one of the presets when applicable
            if ~isempty(presets.(propName))
                if strcmp(type.(propName),'numeric')
                    m = ismember(obj.(propName),cell2mat(presets.(propName)));
                else
                    m = ismember(obj.(propName),presets.(propName));
                end
                if ~m
                    if isa(obj.(propName),'numeric')
                        input_str = num2str(obj.(propName));
                    else
                        input_str = obj.(propName);
                    end
                    if isa(defaults.(propName),'numeric')
                        def_str = num2str(defaults.(propName));
                    else
                        def_str = defaults.(propName);
                    end
                    disp(['Warning: setting ',propName,...
                        ' input value was set to ''',input_str,...
                        ''' which is not one of the acceptable presets. Using ',...
                        def_str,' instead.']);
                    obj.(propName) = defaults.(propName);
                end
            end
        end
        
        % collects list of files to analyze and figures out the 
        % dimensionality of the data to analyze from the first file in the
        % list
        function [obj,fList] = getFListAndNumDim(obj)
        % collects the list of files to analyze
        % figures out the dimensionality of the data to analyze from the first file
        % in the list.

            fList = obj.getFList();
            if isempty(fList)
                obj.fileProcessingMode = 'cancel';
                return
            end
            
            obj.getNumDim(fList);
        end
        
        % figures out the dimensionality of the data to analyze from the first file
        % in flist.
        function obj = getNumDim(obj,fList)
        % figures out the dimensionality of the data to analyze from the first file
        % in fList.
        % enters it as params.numDim

            % quit if no files in the list
            if isempty(fList)
                obj.fileProcessingMode = 'cancel';
                return
            end

            % find the dimensionality of the first image in the list
            [nd,status] = getImageDimensionality(fList{1});
            if strcmp(status,'cancel')
                obj.fileProcessingMode = 'cancel';
                return
            else
                if ismember(obj.fileProcessingMode,...
                        {'singleFileMovie','movieBatch'})
                    % if file is a movie, subtract one to the dimensionality
                    obj.numDim = nd-1;
                else
                    obj.numDim = nd;
                end
            end

            % check that dimensionality is either 2 or 3
            if ~ismember(obj.numDim,[2,3])
                disp(['data is ',num2str(obj.numDim),...
                    'D; should be either 2D or 3D']);
                obj.fileProcessingMode = 'cancel';
                return
            end

        end
        
        % collects list of files to analyze
        function fList = getFList(obj)
            if ismember(obj.fileProcessingMode,{'batch','movieInDir','batchMovie'}) 
                if ~exist(obj.dataFileName,'dir')
                    fList = [];
                    disp(['Could not find data folder ',obj.dataFileName]);
                    return
                else
                    caseSensitive = 0;
                    fList = get_clean_file_list(obj.dataFileName,...
                        obj.inclusionString , obj.exclusionString,...
                        obj.recursive,caseSensitive);
                end
            else
                if ~exist(obj.dataFileName,'file')
                    fList = [];
                    disp(['Could not find data file ',obj.dataFileName]);
                    return
                else
                    fList = {obj.dataFileName};
                end
            end
        end
        
        % converts parameters object into a legacy parameter object
        % (versions < 1.5) or vice versa
        function varargout = legacyParamsConversion(obj,p,alData,direction)
            % direction can be either 'toLegacy' or fromLegacy' to set the
            % conversion direction.
            % p entry is ignored if converting to legacy
            % p entry should be a legacy format object when converting from
            % legacy.
            % alData is ignored if converting from legacy
            % alData should be an airlocalizeData object if converting to
            % legacy.
            % when converting from legacy, the input object is modified (no output);
            % when converting to legacy, the legacy oject is the output.
            v = obj.getVariableCheatSheet;
%             for i=1:numel(v)
%                 disp(v{i});
%             end
            switch direction
                case 'toLegacy'
                    p2 = obj.setLegacyDefaults();
                    f = properties(obj);
                    for i=1:numel(f)
                        if ~strcmp(f{i},'n/a')
                            idx = find(ismember(v(:,1),f{i}));
                            if ~isempty(idx)
                                if numel(idx) == 1
                                    if strcmp(f{i},'threshUnits')
                                        p2.thresh.units = obj.threshUnits;
                                        if strcmp(obj.threshUnits,'SD')
                                            p2.thresh.isLegacySD = 'no';
                                        elseif strcmp(obj.threshUnits,'legacySD')
                                            p2.thresh.isLegacySD = 'yes';
                                        end
                                    elseif strcmp(f{i},'threshLevel')
                                        p2.thresh.level = obj.threshLevel;
                                    elseif strcmp(f{i},'filterLo')
                                        p2.filter.nlo = obj.filterLo;
                                    elseif strcmp(f{i},'filterHi')
                                        p2.filter.nhi = obj.filterHi;
                                    elseif ~strcmp(v{idx,2},'n/a')
                                        p2.(v{idx,2}) = obj.(f{i});
                                    end
                                else
                                    if strcmp(f{i},'psfSigma')
                                        p2.sigma_xy = obj.psfSigma(1);
                                        p2.sigma_z = obj.psfSigma(2);
                                    end 
                                end
                            end
                        end
                    end
                    p2.data_stackname = alData.curFile;
                    varargout{1} = p2;
                case 'fromLegacy'
                    % reset object
                    obj.reset();
                    f = fieldnames(p);
                    for i=1:numel(f)
                        if ~strcmp(f{i},'n/a')
                            idx = find(ismember(v(:,2),f{i}));
                            if ~isempty(idx)
                                if numel(idx) == 1
                                    if strcmp(f{i},'thresh')
                                        obj.threshUnits = p.thresh.units;
                                        obj.threshLevel = p.thresh.level;
                                    elseif strcmp(f{i},'filter')
                                        obj.filterLo = p.filter.nlo;
                                        obj.filterHi = p.filter.nhi;
                                    elseif strcmp(f{i},'sigma_xy')
                                        obj.psfSigma(1) = p.sigma_xy;    
                                    elseif strcmp(f{i},'sigma_z')
                                        obj.psfSigma(2) = p.sigma_z;
                                    elseif ~strcmp(v{idx,1},'n/a')
                                        obj.(v{idx,1}) = p.(f{i});
                                    end
                                end
                            end
                        end
                    end
                otherwise
                    varargout{1} = [];
                    if ~ischar(direction)
                        disp(['Cannot convert params to/from legacy; ',...
                        'input direction was not a character; ',...
                        ' should be set to either to ''toLegacy'' or ''fromLegacy''']);
                        return
                    end
                    disp(['Cannot convert params to/from legacy; ',...
                        'input direction was set to ',direction,...
                        ' should be set to either to ''toLegacy'' or ''fromLegacy''']);
            end
        end

        % returns defaults and accecptable formats for each property
        function [d,t,p,n,section] =  getSettingsDefaultAndPresets(~)
            % outputs defaults and acceptable types/entried for AIRLOCALIZE settings
            % as structures.
            % d: structure holding the default values for each setting used by AIRLOCALIZE
            % t: structure holding the default type for each setting used by AIRLOCALIZE
            % p: structure holding the presets for settings used by AIRLOCALIZE
            % n: structure holding a flag telling whether each setting is absolutely needed for AIRLOCALIZE
                % 0 if parameter is optional (substituted to default without warning)
                % 1 if parameter is always needed
                % 2 if parameter is needed only in some circumstances
                % 3 if parameter is needed but can be substituted to static default with a
                % warning
                % 4 if parameter is needed but can be substituted to dynamic default with a
                % warning

            %% [Files]  

            % path to data file name / folder name for batch analysis
            d.dataFileName = '';
            t.dataFileName = {'char'};
            p.dataFileName = {};
            n.dataFileName = 1;
            section.dataFileName = 'Files';
            
            % path to saving directory
            d.saveDirName = '';
            t.saveDirName = {'char'};
            p.saveDirName = {};
            n.saveDirName = 4;
            section.saveDirName = 'Files';

            % whether to process an individual image or files in batch; whether files are movies or not.    
            d.fileProcessingMode = 'singleFile';
            t.fileProcessingMode = {'char'};
            p.fileProcessingMode = {'singleFile','batch',...
                'movieInDir','singleFileMovie',...
                'batchMovie'};
            n.fileProcessingMode = 1;
            section.fileProcessingMode = 'Files';

            % comma-separated string containing patterns that need to be present in the
            % full file path for it to be analyzed
            d.inclusionString = '';
            t.inclusionString = {'char'};
            p.inclusionString = {};
            n.inclusionString = 0;
            section.inclusionString = 'Files';

            % comma-separated string containing patterns that if present in the full 
            % file path, will exclude the file from analysis
            d.exclusionString = '';
            t.exclusionString = {'char'};
            p.exclusionString = {};
            n.exclusionString = 0;
            section.inclusionString = 'Files';
            
            % set to 1 to explore subdirectories of the parent folder (batch mode
            % only)
            d.recursive = 0;
            t.recursive = {'numeric'};
            p.recursive = {0,1};
            n.recursive = 0;
            section.recursive = 'Files';

            %% [detectionSettings]
            d.psfSigma =  [1,2];
            t.psfSigma =  {'numeric','>0'};
            p.psfSigma =  {};
            n.psfSigma =  1;
            section.psfSigma = 'detectionSettings';

            d.threshUnits = 'absolute' ;
            t.threshUnits = {'char'} ;
            p.threshUnits = {'absolute','SD','legacySD'};
            n.threshUnits = 1 ;
            section.threshUnits = 'detectionSettings';

            d.threshLevel = 1000 ;
            t.threshLevel = {'numeric'};
            p.threshLevel = {} ;
            n.threshLevel = 1 ;
            section.threshLevel = 'detectionSettings';
            
            d.maxSpots = 20000 ;
            t.maxSpots = {'numeric'} ;
            p.maxSpots = {} ;
            n.maxSpots = 3 ;
            section.maxSpots = 'detectionSettings';
            
            d.fitMethod = '3DMaskFull' ;
            t.fitMethod = {'char'} ;
            p.fitMethod = {'3DMaskFull','2DMaskOnLocalMaxProj','3DGaussianFit',...
                '2DGaussianMask','2DGaussianFit'} ;
            n.fitMethod = 4;
            section.fitMethod = 'detectionSettings';
            
            %% [userInput]
            d.setPsfInteractively = 0; 
            t.setPsfInteractively = {'numeric'}; 
            p.setPsfInteractively = {0,1}; 
            n.setPsfInteractively = 0; 
            section.setPsfInteractively = 'userInput';
            
            d.setThreshInteractively = 0; 
            t.setThreshInteractively = {'numeric'}; 
            p.setThreshInteractively = {0,1}; 
            n.setThreshInteractively = 0; 
            section.setThreshInteractively = 'userInput';
            
            d.movieFrameUsedToSetParameters = 1; 
            t.movieFrameUsedToSetParameters = {'numeric','>0','integer'}; 
            p.movieFrameUsedToSetParameters = {}; 
            n.movieFrameUsedToSetParameters = 0; 
            section.movieFrameUsedToSetParameters = 'userInput';
            
            %% [outputImg]
            d.outputSpotsImage = 0;
            t.outputSpotsImage = {'numeric'};
            p.outputSpotsImage = {0,1};
            n.outputSpotsImage = 0;
            section.outputSpotsImage = 'userInput';
            
            d.autoSpotSize = 1;
            t.autoSpotSize = {'numeric'};
            p.autoSpotSize = {0,1};
            n.autoSpotSize = 0;
            section.autoSpotSize = 'userInput';
            
            d.spotSize = [1,2];
            t.spotSize = {'numeric'};
            p.spotSize = {};
            n.spotSize = 2;
            section.spotSize = 'userInput';
            
            d.autoSpotIntensity = 1;
            t.autoSpotIntensity = {'numeric'};
            p.autoSpotIntensity = {0,1};
            n.autoSpotIntensity = 0;
            section.autoSpotIntensity = 'userInput';
            
            d.spotIntensity = 1000;
            t.spotIntensity = {'numeric','>0'};
            p.spotIntensity = {};
            n.spotIntensity = 2;
            section.spotIntensity = 'userInput';
            
            %% [advanced]
            d.minDistBetweenSpots = 2; 
            t.minDistBetweenSpots = {'numeric','>0'}; 
            p.minDistBetweenSpots = {}; 
            n.minDistBetweenSpots = 0;
            section.minDistBetweenSpots = 'advanced';
            
            d.filterLo = 2.15; 
            t.filterLo = {'numeric','>0'}; 
            p.filterLo = {}; 
            n.filterLo = 0;
            section.filterLo = 'advanced';
            
            d.filterHi = 1; 
            t.filterHi = {'numeric','>0'}; 
            p.filterHi = {}; 
            n.filterHi = 0;
            section.filterHi = 'advanced';
            
            d.fittedRegionSize = 3; 
            t.fittedRegionSize = {'numeric','>0'}; 
            p.fittedRegionSize = {}; 
            n.fittedRegionSize = 0; 
            section.fittedRegionSize = 'advanced';
            
            d.tol = 0.01; 
            t.tol = {'numeric','>0'}; 
            p.tol = {}; 
            n.tol = 0; 
            section.tol = 'advanced';
            
            d.bgRegionThickness = 1; 
            t.bgRegionThickness = {'numeric','>0','integer'}; 
            p.bgRegionThickness = {}; 
            n.bgRegionThickness = 0;
            section.bgRegionThickness = 'advanced';
            
            d.psfType = 'integratedGaussian'; 
            t.psfType = {'char'};
            p.psfType = {'integratedGaussianStdZ','integratedGaussian','gaussian'};
            n.psfType = 0;
            section.psfType = 'advanced';
            
            d.maxIterations = 100 ; 
            t.maxIterations = {'numeric','>0'} ; 
            p.maxIterations = 100 ; 
            n.maxIterations = 0 ; 
            section.maxIterations = 'advanced';
            
            d.bgCorrectionMode = 'localPlane';
            t.bgCorrectionMode = {'char'};
            p.bgCorrectionMode = {'localPlane','localMedian','globalMedian'};
            n.bgCorrectionMode = 0;
            section.bgCorrectionMode = 'advanced';
            
            d.useLegacyFormats = 0;
            t.useLegacyFormats = {'numeric'};
            p.useLegacyFormats = {0,1};
            n.useLegacyFormats = 0;
            section.useLegacyFormats = 'advanced';
            
            %% internal use
            d.numDim = 3; % dimensionality is set internally by checking the size of the first image (getNumDim)
            t.numDim = {'numeric'};
            p.numDim = {2,3};
            n.numDim = 4;
            section.numDim = 'internal';
        end
        
        % outputs a structure where each row lists the new property name
        % (v{i,1)) and the old field name (v{i,2}) of params objects for 
        % versions >= 1.6 (new) or <= 1.5 (old).
        function v = getVariableCheatSheet(~)
            v = { ... 
                 'fileProcessingMode','mode';...
                 'dataFileName','data_stackname';...
                 'saveDirName','save_dirname';...
                 'inclusionString','inclusion_string';...
                 'exclusionString','exclusion_string';...
                 'psfSigma' , 'sigma_xy';...
                 'psfSigma' , 'sigma_z';...
                 'threshUnits' , 'thresh.units' ;...
                 'threshLevel' , 'thresh.level';...
                 'maxSpots' , 'max_spots';...
                 'fitMethod', 'fit';...
                 'setPsfInteractively' , 'set_psf_manually';...
                 'setThreshInteractively' , 'set_thresh_manually';...
                 'movieFrameUsedToSetParameters' , 'movie_frame_used_to_set_parameters';...
                 'outputSpotsImage' , 'output_spot_image';...
                 'autoSpotSize' , 'n/a';...
                 'spotSize' , 'n/a';...
                 'autoSpotIntensity' , 'n/a';...
                 'spotIntensity' , 'n/a';...
                 'useDefaults' , 'n/a';...
                 'minDistBetweenSpots' , 'ROIsize';...
                 'filterLo' , 'filter.nlo';...
                 'filterHi' , 'filter.nhi';...
                 'fittedRegionSize', 'cutsize';...
                 'tol' , 'tol';...
                 'bgRegionThickness' , 'thickness';...
                 'psfType' , 'type';...
                 'maxIterations' , 'maxcount';...
                 'bgCorrectionMode' , 'bg_mode';...
                 'n/a', 'cutwidth';...
                 'numdim', 'numDim';};
        end
        
        % outputs a legacy object p with all fields set to defaults
        function p = setLegacyDefaults(~)
            p = struct(...
              'sigma_xy',{},...
                'sigma_z',{},...
                'dx',{},...
                'dy',{},...
                'dz',{},...
                'cutsize',{},...
                'tol',{},...
                'thickness',{},...
                'thresh',{},...
                'max_spots',{},...
                'ROIsize',{},...
                'method',{},...
                'fit',{},...
                'maxcount',{},...
                'data_fname',{},...
                'data_stackname',{},...
                'data_dirname',{},...
                'save_dirname',{});
            
                p(1).numdim = 3;

                % set default values
                p.sigma_xy = 2;
                p.sigma_z = 2;
                p.dx = 64;
                p.dy = 64;
                p.dz = 200;
                p.cutsize = 3;
                p.tol = 0.01;
                p.thickness = 1;
                p.thresh(1).level = 3;
                p.thresh.units = 'SD';
                p.thresh.sd = 0;
                p.filter.nlo = 2.15;
                p.filter.nhi = 1;
                p.max_spots = 200000; 
                p.ROIsize = 2;
                p.fit = '3d mask';
                p.maxcount = 100;
                p.data_fname = [];
                p.data_stackname = [];
                p.data_dirname = [];
                p.save_dirname = [];
                p.type = 'integrated gaussian';
                p.set_psf_manually = 1;
                p.set_thresh_manually = 1;
                p.data = 0;
                p.smooth = 0;
                p.bg_mode = 'local plane';
                p.mode = 'null';
        end
    end
end

