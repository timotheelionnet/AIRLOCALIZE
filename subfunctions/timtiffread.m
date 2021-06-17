function [stack,nSlices] = timtiffread(varargin)
if nargin == 0 
    disp('no arguments loaded into tiff reader, cannot read');
    stack = [];
    nSlices = 0;
    return
    
else 
    try
        stack = tiffread5(varargin{1}); % this sometimes has issues reading files
        nSlices = 1;
        if ndims(stack) == 3
            nSlices = size(stack,3);
        end
        disp(' loaded image file with tifread');
    catch
        [stack,nSlices] = readTifStackWithImRead(varargin{1}); % more stable reader but sometimes slower
        disp(' loaded image file with imRead');
    end


end

