function [params,alData] = set_psf_size_manually(params,alData)
    
    % retrieve current image into alData object
    overwrite = 0;
    alData.retrieveImg(overwrite);
    if ismember(params.fileProcessingMode,...
            {'singleFileMovie','batchMovie'})
        %alData.img = ...
        %    alData.img(:,:,params.movieFrameUsedToSetParameters);
    end
    
    res = explore_stack_v6(params,alData);
    
    if max(res(:)) == 0
        params.fileProcessingMode = 'cancel';
        disp('you did not choose any spot.');
        return
    end
    
    if params.numDim == 3
        params.psfSigma(1) = mean(res(:,6));
        params.psfSigma(2) = mean(res(:,7));
    elseif params.numDim == 2    
        params.psfSigma = mean(res(:,5));
    end

end
