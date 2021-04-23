function flist = get_clean_file_list(inputdir,inc, exc,recursive,case_sensitive)
%input: a directory path (inputdir)
%a cell array containing strings of characters that are included in the
%files to be kept (inc) - set to empty if unnecessary
%only file names containing ALL strings from the cell array will be retained

%a cell array containing strings of characters that are included in the
%files to be excluded (exc)- set to empty if unnecessary
%all file names containing at least ONE string from the cell array will be discarded

%a number (0 or 1) that sets the recursive option (1 = recursive; 0 = root
%folder only)
%case sensitive: set to 0 or 1
%outputs a list of all the file names (full file name including path)

if recursive == 0
   flist = dir( inputdir);
   %remove directories
   flist = flist(~[flist.isdir]);
   dircell = cell(size(flist));
   dircell(:) = {inputdir};
   fnames = struct2cell(flist);
   fnames = fnames(1,:);
   fnames = fnames';
   flist = cellfun(@fullfile,dircell,fnames,'UniformOutput',0);
else
    %get all files
   flist = dirrec(inputdir);
end

if ~isempty(inc)
    for i=1:numel(inc)
        tmpinc = inc{i};
        if ~isempty(inc{i})
            if ~case_sensitive
                x = regexpi(flist,tmpinc);
            else
                x = regexp(flist,tmpinc);
            end
            idx = ~(cellfun(@isempty,x));
            flist = flist(logical(idx));
        end
    end  
end

if ~isempty(exc)
    for i=1:numel(exc)
        tmpexc = exc{i};
        if ~isempty(exc{i})
            if ~case_sensitive
                x = regexpi(flist,tmpexc);
            else
                x = regexp(flist,tmpexc);
            end
            idx = cellfun(@isempty,x);
            flist = flist(logical(idx));
        end
    end  
end
if isrow(flist)
    flist = flist';
end
end