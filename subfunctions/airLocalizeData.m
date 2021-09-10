classdef airLocalizeData < handle
    
    % object that manages access to data for Airlocalize
    
    properties (GetAccess = 'public', SetAccess = 'public')
        fList; % holds the file list
        curFile; % holds the name of the current file
        fileIdx; % holds the index of the current file in fList
        img; % holds the raw image/stack data of the current file
        smooth; % holds the smoothed image/stack data of the current file
        isMovie; % 1 if files in the list are movies, 0 if files are individual images
        nFrames; % number of frames in current file (if movie)
        curFrame; % current frame in current file (if movie)
    end
    
    methods
        
        % constructor
        function obj = airLocalizeData()
             obj.reset();
        end
        
        function obj = reset(obj)
            obj.fList = {}; 
            obj.curFile = ''; 
            obj.fileIdx = 0; 
            obj.img = []; 
            obj.smooth = [];
            obj.isMovie = 0;
            obj.nFrames = 0;
            obj.curFrame = 0;
            end
        
        function obj = resetCurrentFile(obj)
            % resets everything but the file list and the isMovie flag.
            obj.curFile = ''; 
            obj.fileIdx = 0; 
            obj.img = []; 
            obj.smooth = [];
            obj.nFrames = 0;
            obj.curFrame = 0;
        end
        
        % generates a list of files to analyze
        % basd on dataFileName and fileProcessingMode properties of params
        % airlocalizeParams object.
        function obj = setFListFromParams(obj,params)
            obj.reset();
            if ismember(params.fileProcessingMode,...
                    {'batch','movieInDir','batchMovie'})
                if ~exist(params.dataFileName,'dir')
                    obj.fList = obj.setFList({});
                    disp(['Could not find data folder ',params.dataFileName]);
                    return
                else
                    caseSensitive = 0;
                    
                    %format comma-separated list into cell arrays of substrings
                    inclusionString = parse_and_format_string_list(params.inclusionString);
                    exclusionString = parse_and_format_string_list(params.exclusionString);
                    imgExt = getAllowedImageFormatExtensions();
                    fl = {};
                    for i=1:numel(imgExt)
                        tmpfl = get_clean_file_list(params.dataFileName,...
                            [inclusionString,imgExt(i)] , exclusionString,...
                            params.recursive,caseSensitive);
                        idx = ones(numel(tmpfl),1);
                        for j=1:numel(tmpfl) %making sure that the files have a legit image extension
                            [~,~,e] = fileparts(tmpfl{j});
                            if ~ismember(e,imgExt)
                                idx(j) = 0;
                            end
                        end
                        tmpfl = tmpfl(logical(idx));
                        fl = [fl;tmpfl];
                    end
                    fl = unique(fl);
                end
            else
                if ~exist(params.dataFileName,'file')
                    obj.fList = obj.setFList([]);
                    disp(['Could not find data file ',params.dataFileName]);
                    return
                else
                    fl = {params.dataFileName};
                end
            end
            obj.setFList(fl);
            if ismember(params.fileProcessingMode,...
                    {'singleMovie','batchMovie'})
                obj.isMovie = 1;
            end
        end
        
        % sets a file list as the fList property
        function obj = setFList(obj,fList)
            obj.reset();
            if iscell(fList)
                obj.fList = fList;
            else 
                disp(['input file list has wrong format;',...
                    ' resetting airlocalize data object']);
            end
        end
        
        % outputs the file list
        function fList = getFList(obj)
            fList = obj.fList;
        end
        
        % sets the index of the current file if compatible with fList size
        function obj = setFileIdx(obj,idx)
            if idx > numel(obj.fList) || idx <= 0
                disp(['cannot access file index ',num2str(idx),...
                    '; file list has ',num2str(numel(obj.fList)),...
                    ' entries']);
                return
            end
            if idx ~= obj.fileIdx
                obj.resetCurrentFile;
                obj.fileIdx = idx;
                obj.curFile = obj.fList{idx};
            else
                if ~strcmp(obj.curFile,obj.fList{idx})
                    obj.resetCurrentFile;
                    obj.fileIdx = idx;
                    obj.curFile = obj.fList{idx};
                end
            end
        end
        
        function idx = getFileIdx(obj)
            idx = obj.fileIdx;
        end
        
        function obj = setCurFile(obj,fileName)
            if ~ismember(fileName,obj.fList)
                disp(['desired file ',fileName,...
                    ' is not part of existing file list; ',...
                    'cannot set as current.']);
            else
                idx = find(ismember(obj.fList,fileName));
                obj.setFileIdx(idx);
            end
        end
        
        function fileName = getCurFile(obj)
            fileName = obj.curFile;
        end
        
        function obj = setNFrames(obj)
            overwrite = 0;
            if ~obj.isFileIndexImgLoaded(obj.fileIdx)
                obj.retrieveImg(overwrite);
            else
                s = size(obj.img);
                if numel(s) == 2
                    obj.nFrames = 1;
                else
                    obj.nFrames = s(3);
                end
            end
        end
        
        function status = setCurFrame(obj,newFrame)
            if newFrame == obj.curFrame
                status = 1;
                return
            end
            if round(newFrame) ~= newFrame  || newFrame < 1 || ...
                    newFrame > obj.nFrames
                status = 0;
            else
                obj.curFrame = newFrame;
            end
        end
        
        function isLoaded = isFileIndexImgLoaded(obj,idx)
            isLoaded = 0;
            if obj.fileIdx == idx && idx > 0 && idx <= numel(obj.fList) ...
                && ~isempty(obj.curFile) && ~isempty(obj.img)
                if numel(obj.fList) >= idx
                    if strcmp(obj.curFile,obj.fList{idx})
                        isLoaded = 1;
                    end
                end
            end    
        end
        
        function isLoaded = isFileNameImgLoaded(obj,fileName)
            isLoaded = 0;
            if obj.fileIdx > 0 && obj.fileIdx <= numel(obj.fList)
                if strcmp(obj.fList{obj.fileIdx},fileName) ...
                    && ~isempty(obj.curFile) && ~isempty(obj.img)
                    if strcmp(obj.curFile,obj.fList{idx})
                        isLoaded = 1;
                    end
                end
            end    
        end
        
        function isLoaded = isFileIndexSmoothedImgLoaded(obj,idx)
            isLoaded = 0;
            if obj.fileIdx == idx && idx > 0 && idx <= numel(obj.fList) ...
                && ~isempty(obj.curFile) && ~isempty(obj.smooth)
                if numel(obj.fList) >= idx
                    if strcmp(obj.curFile,obj.fList{idx})
                        isLoaded = 1;
                    end
                end
            end    
        end
        
        function isLoaded = isFileNameSmoothLoaded(obj,fileName)
            isLoaded = 0;
            if obj.fileIdx > 0 && obj.fileIdx <= numel(obj.fList)
                if strcmp(obj.fList{obj.fileIdx},fileName) ...
                    && ~isempty(obj.curFile) && ~isempty(obj.smooth)
                    if strcmp(obj.curFile,obj.fList{idx})
                        isLoaded = 1;
                    end
                end
            end    
        end
        
        function obj = retrieveImg(obj,overwrite)
            if obj.fileIdx < 0 || obj.fileIdx > numel(obj.fList)
                obj.img = [];
                return
            else
                if ~obj.isFileIndexImgLoaded(obj.fileIdx) || overwrite
                    [obj.img,nSlices] = ...
                        timtiffread(obj.fList{obj.fileIdx});
                    obj.smooth = [];
                    obj.curFrame = 0;
                    if obj.isMovie
                        obj.nFrames = nSlices;
                    else
                        obj.nFrames = 1;
                    end
                end
            end
        end
        
        function obj = retrieveSmoothedImg(obj,params,overwrite)
            if obj.fileIdx < 0 || obj.fileIdx > numel(obj.fList)
                obj.smooth = [];
                return
            else
                if ~obj.isFileIndexSmoothedImgLoaded(obj.fileIdx) || overwrite
                    if ~obj.isFileIndexImgLoaded(obj.fileIdx) || overwrite
                        obj.retrieveImg(overwrite);
                    end
                    obj.smooth = smooth_image_and_subtract_background7(...
                        obj.img,params);
                end
            end
        end
    end
end
