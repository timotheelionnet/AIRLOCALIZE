%%
% root name for AIRLOCALIZE repo
rootFolderShortName = 'AIRLOCALIZE';

% where the subfunctions are
subfunctionsFolderName = 'subfunctions';

% where the obsolete subfunctions will be dumped
obsoleteFolderName = 'obsolete';

%% move to AIRLOCALIZE root
cd('..');
curFolder = pwd;
rootString = strsplit(pwd,filesep);
rootString = rootString{end};
if ~strcmp(rootString,rootFolderShortName)
    disp('You should set your working folder to ',rootFolderShortName,...
        '/management before running this code.');
    return
end

%% prep the obsolete repository

if ~exist(obsoleteFolderName,'dir')
    mkdir(obsoleteFolderName);
end

%% collect list of needed subfunctions
[neededFList,~] = matlab.codetools.requiredFilesAndProducts('AIRLOCALIZE.m');
neededFList = neededFList';

%% collect current subfunction list
currentFList = get_clean_file_list(subfunctionsFolderName,...
            {''}, {''},...
            1,1);
        
% concatenate with the root folder to obtain the full path name
for i=1:numel(currentFList)
    currentFList{i} = [ curFolder , filesep, currentFList{i}];
end

%% collect list of obsolete functions and move them
obsoleteFList = setdiff(currentFList,neededFList);
disp('List of files moved to obsolete folder:');
obsoleteFList
for i=1:numel(obsoleteFList)
    [~,f,e] = fileparts(obsoleteFList{i});
    destName = fullfile(obsoleteFolderName,[f,e]);
    movefile( obsoleteFList{i}, destName);
end



        