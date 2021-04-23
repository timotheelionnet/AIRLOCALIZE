function [pl,stack_bg_corr,new_ctr,ROIlimits] = subtract_median2(...
    alData,spotCenter,cutWidth,thickness,ROIsize,median_type)

if ~alData.isMovie
    imSize = size(alData.img);
else
    imSize = size(alData.img(:,:,alData.curFrame));
end

[new_ctr,ROIlimits] = compute_ROI_boundaries2(...
    imSize,spotCenter,cutWidth,thickness,ROIsize);

if ndims(alData.img) == 3
    stack_bg_corr = alData.img(...
        ROIlimits(1,1):ROIlimits(2,1),...
        ROIlimits(1,2):ROIlimits(2,2),...
        ROIlimits(1,3):ROIlimits(2,3));
else
    if ~alData.isMovie
        stack_bg_corr = alData.img(...
            ROIlimits(1,1):ROIlimits(2,1),...
            ROIlimits(1,2):ROIlimits(2,2));
    else
        stack_bg_corr = alData.img(...
            ROIlimits(1,1):ROIlimits(2,1),...
            ROIlimits(1,2):ROIlimits(2,2),...
            alData.curFrame);
    end
end

if strcmp(median_type,'local')
    stack_bg_corr = stack_bg_corr - median(double(stack_bg_corr(:)));
elseif strcmp(median_type,'global')
    if ~alData.isMovie
        stack_bg_corr = stack_bg_corr - median(double(alData.img(:)));
    else
        stack_bg_corr = stack_bg_corr ...
            - median(double(alData.img(:,:,alData.curFrame),'all'));
    end
end

pl = [];
end