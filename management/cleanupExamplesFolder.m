function cleanupExamplesFolder()
    addpath('../subfunctions');
    exFolderName = '../examples';
    recursive = 1;
    caseSensitive = 0;
    keyWordsToDelete = {'_spots','.loc3','.par','.loc4','.det'};
    keyWordsToKeep = {''};
    % loop through key words to delete and build delete list
    listofFilesToClean = {};
    for i=1:numel(keyWordsToDelete)
        fl = get_clean_file_list(exFolderName,...
            keyWordsToDelete(i), keyWordsToKeep,...
            recursive,caseSensitive);
        listofFilesToClean = [listofFilesToClean;fl];
    end
    % remove duplicates
    listofFilesToClean = unique(listofFilesToClean)
    
    % delete selected files
    for i=1:numel(listofFilesToClean)
        delete(listofFilesToClean{i});
    end
    disp(['deleted ',num2str(numel(listofFilesToClean)),' files.']);
end