function save_as_tiff(varargin)
%normalize = 1 sets the max intensity to be 65536
%otherwise converting the image to int16 format without normalization
%save_as_tiff(img,filename)
%save_as_tiff(img,filename,normalize)
if nargin <2
    disp('I need at least an image and a file name');
    return;
elseif nargin == 2
    img = varargin{1};
    filename = varargin{2};
    normalize = 0;
elseif nargin == 3
    img = varargin{1};
    filename = varargin{2};
    normalize = varargin{3};
elseif nargin > 3
    disp('warning: too many input arguments');
    img = varargin{1};
    filename = varargin{2};
    normalize = varargin{3};
end


if ndims(img) == 2 || (ndims(img) == 3 && size(img,3) == 1)
    %% 2d img
    if(normalize==1)
        img = (img-min(min(img)))/max(max(img))*65536;
    end
    img = uint16(img);

    imwrite(img,filename,'Compression','none');
    return
    
elseif (ndims(img) == 3 && size(img,3) ~= 1)
    %% 3d stack 
    if(normalize==1)
        img = (img-min(img(:)))/(max(img(:)) - min(img(:)) )*65536;
    end
    
    img = uint16(img);
    
    [nx ny nz] = size(img);
    imwrite(img(1:nx,1:ny,1),filename,'Compression','none','WriteMode','overwrite' );
    for i = 2:nz
        imwrite(img(1:nx,1:ny,i),filename,'Compression','none','WriteMode','append');
    end
    return
    
end

end