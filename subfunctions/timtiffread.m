function [stack,nSlices] = timtiffread(varargin)

[stack,nSlices] = readTifStackWithImRead(varargin{1});

end

