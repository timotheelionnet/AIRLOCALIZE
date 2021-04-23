function maxima = find_isolated_maxima_clean_Img(img,params)
% finds local maxima with positive intensity in the input image img

% params is an airlocalizeParams object that needs to contain the following fields: 
    % minDistBetweenSpots: minimum allowed distance between local maxima

% outputs maxima: a list of local maxima pixel coordinates
% [ x y I ] in 2D
% [ x y z I ] in 3D

% generate dummy airlocalizeData object that will hold the image
alData = airLocalizeData();
alData.smooth = img;
alData.img = img;
alData.isMovie = 0;

% generate dummy airlocalizeParams object that will hold the parameters
p = airLocalizeParams;
p.threshLevel = 0;
p.threshUnits = 'absolute';
p.minDistBetweenSpots = params.minDistBetweenSpots;
p.numDim = params.numDim;
p.maxSpots = params.maxSpots;

verbose = 0;
maxima = find_isolated_maxima_clean3(alData,p,verbose);    
end