function res = is_recognized_img_format(fname)
% this functions checks that a file has the expected format

    if ~ischar(fname)
        res = 0;
        return;
    else
        [~,~,ext] = fileparts(fname);
        
        %%modify the next line to add a recognized format
        res = strcmp(ext,'.tiff') || strcmp(ext,'.TIFF') || strcmp(ext,'.tif') || strcmp(ext,'.stk') || strcmp(ext,'.TIF') || strcmp(ext,'.STK');
    end
    
end