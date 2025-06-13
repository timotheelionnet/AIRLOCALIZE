%% example paths
fl1 = get_clean_file_list('/Volumes/lionnt01lab/lionnt01labspace/Minghan Yang/Data/CX7/Bursting/20250324-p300-HIV-Triptolide/AcquireOnly.V3_03-24-25_04;05;12',{'.TIFF'} , {[]},1,0);
fl2 = get_clean_file_list('/Volumes/lionnt01lab/lionnt01labspace/Minghan Yang/Data/CX7/Bursting/20250324-p300-HIV-Triptolide/mask',{'.tif'} , {[]},1,0);

%%
tic;
fileNamesOnly = 0;
fullList = [fl1;fl2];
shortList = trimPaths(fullList); % eliminate common root dir before trying to match files
% [matchedPairs, similarityScores] = match_file_pairs(...
%     shortList(1:numel(fl1)), shortList(numel(fl1)+1:end),fileNamesOnly);
[matchedPairs, similarityScores] = match_file_pairs(fl1, fl2,fileNamesOnly);

pairList = [fl1(matchedPairs(:,1)), fl2(matchedPairs(:,2)) ];
toc

%%
function shortPaths = trimPaths(filePaths)
    % trim common prefix and suffix from file list
    shortPaths = eliminateCommonPrefix(filePaths);
    shortPaths = eliminateCommonSuffix(shortPaths);
end

%%
function shortPaths = eliminateCommonPrefix(filePaths)
% from a list of file path filePaths, eliminate any prefix that is common
% to all files.
    % Convert to char array for line-wise comparison
    paths = char(filePaths);  % n-by-m character array
    minPath = min(paths);
    maxPath = max(paths);
    
    % Compare character by character
    commonLength = find(minPath ~= maxPath, 1, 'first') - 1;
    if isempty(commonLength)
        commonLength = length(minPath);  % Full match
    end
    
    commonPrefix = minPath(1:commonLength);
    
    shortPaths = cellfun(@(p) p(commonLength+1:end), filePaths, 'UniformOutput', false);

end

%%
function shortPaths = eliminateCommonSuffix(filePaths)
% from a list of file path filePaths, eliminate any suffix that is common
% to all files.
    reversedPaths = cellfun(@(p) fliplr(p), filePaths, 'UniformOutput', false);
    shortRevPaths = eliminateCommonPrefix(reversedPaths);
    shortPaths = cellfun(@(p) fliplr(p), shortRevPaths, 'UniformOutput', false);
end

function [matchedPairs, similarityScores] = match_file_pairs(fileList1, fileList2,fileNamesOnly)
% MATCH_FILE_PAIRS Matches files from two lists based on filename similarity.
%   [matchedPairs, similarityScores] = match_file_pairs(fileList1, fileList2)
%   returns a list of index pairs and similarity scores based on token overlap.
%
%   Inputs:
%     - fileList1, fileList2: cell arrays of file paths
%      - set fileNamesOnly to 1 to only perform the matching on the filenames (ignoring the directory)
%
%   Outputs:
%     - matchedPairs: n-by-2 array, rows are [i j] index pairs from list1 and list2
%     - similarityScores: similarity score for each matched pair
    
    if fileNamesOnly
        % Extract filenames
        [~, names1, ext1] = cellfun(@fileparts, fileList1, 'UniformOutput', false);
        [~, names2, ext2] = cellfun(@fileparts, fileList2, 'UniformOutput', false);
        names1 = strcat(names1, ext1);
        names2 = strcat(names2, ext2);
    else
        % keep full path
        names1 = fileList1;
        names2 = fileList2;
    end

    % Initialize similarity matrix
    N1 = numel(names1);
    N2 = numel(names2);
    similarityMatrix = zeros(N1, N2);

    % Compute similarity as number of shared tokens
    for i = 1:N1
        tokens1 = split(lower(names1{i}), {'_', '-', '.', ' '});
        for j = 1:N2
            tokens2 = split(lower(names2{j}), {'_', '-', '.', ' '});
            shared = intersect(tokens1, tokens2);
            similarityMatrix(i, j) = numel(shared);
        end
    end

    % Convert to cost matrix (negative for matchpairs)
    costMatrix = -similarityMatrix;

    % Match using Hungarian algorithm
    [matchedPairs, costs] = matchpairs(costMatrix, -1);

    % Recover similarity scores (invert cost)
    similarityScores = similarityMatrix(sub2ind(size(similarityMatrix), ...
        matchedPairs(:,1), matchedPairs(:,2)));
end
