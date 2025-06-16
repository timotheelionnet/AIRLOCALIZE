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