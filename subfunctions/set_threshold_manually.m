function params = set_threshold_manually(params,alData)

    % retrieve current smoothed image into alData object
    overwrite = 0;
    alData.retrieveSmoothedImg(params,overwrite);
    
    params = explore_stack_2channels4(params,alData); 
    
end