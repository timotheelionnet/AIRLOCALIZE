function [loc,nWrong,nDouble,dummyImg] = clean_up_spots_array_clean4(...
    loc,imSize,params)

%series of checks on the detected spots to remove aberrant detections
if params.numDim == 3
    nx = imSize(1);
    ny = imSize(2);
    nz = imSize(3);
    Icolnum = 4;
elseif params.numDim == 2
    nx = imSize(1);
    ny = imSize(2);
    Icolnum = 3;
end

%% removing the rows corresponding to iterations of the gaussian mask that did not converge 
%(signaled in the gaussian mask algorithm by [x y z N] = [-1 -1 -1 -1])

not_wrong = logical((loc(:,Icolnum) ~= -1)...
    .*(loc(:,1) >= 0).*(loc(:,1) < nx)...
    .*(loc(:,2) >= 0 ).*(loc(:,2) < ny ));

if params.numDim == 3
    not_wrong = logical(not_wrong.*( loc(:,3) >= 0 ).*( loc(:,3) < nz));
end

nWrong = size(loc,1) - sum(not_wrong);
loc = loc(not_wrong,:);
clear('not_wrong');

%% removing residual negative intensity spots
loc = loc(loc(:,Icolnum)>0,:);

%% removing duplicate spots that converged within n_exclusion pixel from each other

%I create a stack/img dummyImg in which each spot is sent to its pixel, with its value
%corresponding to its intensity rank (1= least intense spot)
loc = sortrows(loc,Icolnum,'ascend');
dummyImg = zeros(imSize);
if params.numDim == 3
    linIdx = sub2ind(imSize,ceil(loc(:,1)),ceil(loc(:,2)),ceil(loc(:,3)));
else
    linIdx = sub2ind(imSize,ceil(loc(:,1)),ceil(loc(:,2)));
end

% flipping the order of the list to ensure the brightest spots overwrites
% any dimmer one(s) that might fall onto the same pixel
dummyImg(linIdx) = 1:size(loc,1); 

%I check whether each spot is the maximal ranking within an ROI of
%size n_exclusion, otherwise I discard it
maxima = find_isolated_maxima_clean_Img(dummyImg,params);
% flip the last column of maxima to recover ordered indices
%maxima(:,end) = size(loc,1)-maxima(:,end)+1;
nDouble = size(loc,1)-size(maxima,1);

loc = loc(maxima(:,end),:);
loc = sortrows(loc,Icolnum,'descend');

end
