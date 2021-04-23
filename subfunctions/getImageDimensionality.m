function [numDim,status] = getImageDimensionality(imgName)
    numDim = 0;
    if ~exist(imgName,'file') 
        disp(['cannot find image file ',imgName]);
        disp('cannot get dimensionality');
        status = 'cancel';
        return
    end
    
    if ~is_recognized_img_format(imgName)
        disp(['image file ',imgName,' is not a recognized format;']);
        disp('cannot get dimensionality');
        status = 'cancel';
        return
    end
    
    tifInfo = imfinfo(imgName);
    if size(tifInfo, 1) == 1
        numDim = 2;
        status = 'success';
    elseif size(tifInfo, 1) > 1
        numDim = 3;
        status = 'success';
    else
        disp(['image file ',imgName,' is not a recognized format;']);
        disp('cannot get dimensionality');
        status = 'cancel';
        return
    end
end