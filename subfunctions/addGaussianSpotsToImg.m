function img = addGaussianSpotsToImg(imgSize,loc,spotSize)
% imgSize: size of the output image either [nx ny] or [ny ny nz] / [nx ny
% nframes]
% loc: loc array wit following columns: 

numDim = numel(imgSize);
img = zeros(imgSize);

for i=1:size(loc,1)
    
    % bounding box coordinates around spot
    [bmin,bmax] = boundingBoxCoordinates(...
        imgSize,loc(i,1:numDim+1),spotSize);

     % calculate gaussian intensity across box & add to img 1
    img = addGaussianIntensityToImg(...
        img,bmin,bmax,loc(i,1:numDim),loc(i,numDim+1),spotSize);
end

end

function img = addGaussianIntensityToImg(img,bmin,bmax,spotPosition,...
    spotIntensity,spotSize)
    % add offset to recenter Airlocalize coordinates onto pixels
    spotPosition(1) = spotPosition(1)+0.5;
    spotPosition(2) = spotPosition(2)+0.5;
    
    if ndims(img) == 2
        % create meshgrid swapping x/y to match image coordinates
        [yy,xx] = ...
            meshgrid(bmin(2):bmax(2),bmin(1):bmax(1));
        
        % compute gaussian intensity
        g = exp(-(xx - spotPosition(1)).^2/(2*spotSize(1)^2));
        g = g .* ...
            exp(-(yy - spotPosition(2)).^2/(2*spotSize(1)^2));
        g = g * spotIntensity / sum(g,'all');
        
        img(bmin(1):bmax(1),bmin(2):bmax(2)) = ...
            img(bmin(1):bmax(1),bmin(2):bmax(2)) + g;
    else
        % create meshgrid swapping x/y to match image coordinates
        [yy,xx,zz] = ...
            meshgrid(bmin(2):bmax(2),bmin(1):bmax(1),bmin(3):bmax(3));
        
        % compute gaussian intensity
        g = exp(-(xx - spotPosition(1)).^2/(2*spotSize(1)^2));
        g = g .* ...
            exp(-(yy - spotPosition(2)).^2/(2*spotSize(1)^2));
        g = g .* ...
            exp(-(zz - spotPosition(3)).^2/(2*spotSize(2)^2));
        g = g * spotIntensity / sum(g,'all');
        
        img(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3)) = ...
            img(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3)) + g;
    end
end

function [bmin,bmax] = boundingBoxCoordinates(...
    imgSize,spotPosition,spotSize)

    % coordinates of bounding box for spots in pixel units
    spotSize = 3* spotSize; % taking 3 times the sigma of the gaussian

    bmin(1) = ceil(spotPosition(1) - spotSize(1));
    bmin(1) = max( min(bmin(1),imgSize(1)), 1);
    bmin(2) = ceil(spotPosition(2) - spotSize(1));
    bmin(2) = max( min(bmin(2),imgSize(2)), 1);
    bmax(1) = ceil(spotPosition(1) + spotSize(1));
    bmax(1) = max( min(bmax(1),imgSize(1)), 1);
    bmax(2) = ceil(spotPosition(2) + spotSize(1));
    bmax(2) = max( min(bmax(2),imgSize(2)), 1);
    if numel(imgSize) == 3
        bmin(3) = ceil(spotPosition(3) - spotSize(2));
        bmin(3) = max( min(bmin(3),imgSize(3)), 1);
        bmax(3) = ceil(spotPosition(3) + spotSize(2));
        bmax(3) = max( min(bmax(3),imgSize(3)), 1);
    end
end
    