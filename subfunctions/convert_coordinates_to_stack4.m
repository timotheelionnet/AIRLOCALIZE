function stack = convert_coordinates_to_stack4(varargin)
% convert_coordinates_to_stack4(loc,imSize)
% convert_coordinates_to_stack4(loc,imSize,'ValMode',<valmode>)
% convert_coordinates_to_stack4(loc,imSize,'ValMode',<valmode>,...
%   'Size',<size>)
% convert_coordinates_to_stack4(loc,imSize,'ValMode',<valmode>,...
%   'Size',<size>,'Additive',<additive>)

% valmode optional argument
% valmode = 'value' (default): build a stack with the values corresponding to the last
% column of the array [x y z I] or [x y I]

% valmode = 'ones': positions indicated by the input array are set to ones
%   [x y z I] or [x y I]

% size : size of the square around the center of the spot 
    % I take the closest odd integer from that number

% Additive =1 : adds values of overlapping spots (default)
% =0: overwrites in case of overlap

%%
if nargin < 2
    stack = [];
    return;
end

%% parsing arguments for loc coordinates and the dimension of the stack
loc = varargin{1};
imSize = varargin{2};
numDim = numel(imSize);

if ~(numDim == 2 || numDim == 3)
    stack = [];
    return
end

if size(loc,2) < numDim
    stack = [];
    return
end

%% parsing argument for options
i = 1;
if numDim == 3
    spotSize = [1,1];
else
    spotSize = 1;
end
ValMode = 'ones';
additive = 1;
while i<= nargin
    if isa(varargin{i},'char')
        
        if strcmp(varargin{i},'ValMode')
            if i < nargin
                ValMode = varargin{i+1};
                i=i+2;
            end
            
        elseif strcmp(varargin{i},'Size')
            if i < nargin
                spotSize = varargin{i+1};
                spotSize = round(spotSize);
                i=i+2;
            end
        elseif strcmp(varargin{i},'Additive')
            if i < nargin
                additive = varargin{i+1};
                i=i+2;
            end    
        else
            i = i+1;
        end
    else 
        i = i+1;
    end
end

if ~isa(spotSize,'numeric')
    spotSize = 1;
else
    % id Size is even, take the next odd number 
    % (so that the square is centered on the position
    % of the spot
    if mod(spotSize,2)==0
        spotSize = spotSize+1;
    end
end

if ~isa(additive,'numeric') && ~isa(additive,'logical')
    additive = 1;
else
    if additive ~= 0 && additive ~= 1
        additive = 1;
    end
end

if ~isa(ValMode,'char')
    ValMode = 'ones';
end

%% generating the stack / image
stack = fill_stack(loc,imSize,numDim,ValMode,spotSize,additive);

end

function stack = fill_stack(loc,imSize,numDim,ValMode,spotSize,additive)
    
    stack = zeros(imSize);
    
    % assign each spot to its integer pixel coordinate
    loc(:,1:2) = ceil(loc(:,1:2));
    if numDim== 3
        loc(:,3) = round(loc(:,3));
    end
    
    % eliminate spots that fall out of the image
    outOfImg = (loc(:,1) > imSize(1)) | (loc(:,1) < 1);
    outOfImg = outOfImg | (loc(:,2) > imSize(2)) | (loc(:,2) < 1);
    if numDim == 3
        outOfImg = outOfImg | (loc(:,3) > imSize(3)) | (loc(:,3) < 1);
    end
    
    loc = loc(~outOfImg,:);
    
    for i = 1:size(loc,1)
        if  strcmp(ValMode,'value') && size(loc,2) >= numDim +1
            iValue = loc(i,numDim+1);
        else
            iValue = 1;
        end
        stack = add_spot_to_img(stack,loc(i,1:numDim),iValue,...
            spotSize,additive);
    end
    
    if strcmp(ValMode,'ones')  
        stack = 255*stack;
    end
end

function img = add_spot_to_img(img,centerPos,iValue,spotSize,additive)

imSize = size(img);

dx = ceil(spotSize(1)/2);
x1 = max(1,centerPos(1) - dx);
x2 = min(imSize(1),centerPos(1) + dx);
y1 = max(1,centerPos(2) - dx);
y2 = min(imSize(2),centerPos(2) + dx);
if ndims(img) == 3
    dz = floor(spotSize(2)/2);
    z1 = max(1,centerPos(3) - dz);
    z2 = min(imSize(3),centerPos(3) + dz);
    if additive == 1
        img(x1:x2,y1:y2,z1:z2) = img(x1:x2,y1:y2,z1:z2)+ iValue;
    else
        img(x1:x2,y1:y2,z1:z2) = iValue;
    end
else
    if additive == 1
        img(x1:x2,y1:y2) = img(x1:x2,y1:y2)+ iValue;
    else
        img(x1:x2,y1:y2) = iValue;
    end
end 


end
