function spotsImg = generateSpotsImage(params,loc,locVars,alData)
    
    % if not set in params, set the saving format to the default
    if ~isprop(params,'useLegacyFormats')
        defaults = params.getSettingsDefaultAndPresets();
        params.useLegacyFormats = defaults.useLegacyFormats;
    end
    
    % set spot size to automatic value if needed
    if params.autoSpotSize
        params.spotSize = params.psfSigma;
    end

    if params.useLegacyFormats
        % if legacy format, generate cubes (squares in 2D) around each spot   
        if ~alData.isMovie
            spotsImg = convert_coordinates_to_stack4(loc,...
                size(alData.img),'ValMode','ones',...
                'Size',params.spotSize,'Additive',0);
        else
            spotsImg = zeros(size(alData.img));
            frameSize = size(alData.img);
            frameSize = frameSize(1:2);
            frameColNum = ismember(locVars,'frame');
            for i=1:alData.nFrames
                tmpLoc = loc(loc(:,frameColNum)==i,:);
                spotsImg(:,:,i) = convert_coordinates_to_stack4(tmpLoc,...
                frameSize,'ValMode','ones',...
                'Size',params.spotSize,'Additive',0);
            end
        end
    else
        % if non-legacy format, generate a gaussian distribution around each
        % spot
        if ~params.autoSpotIntensity
            if params.spotIntensity > 0
                loc(:,params.numDim+1) = params.spotIntensity;
            end
        end
        if ~alData.isMovie
            spotsImg = addGaussianSpotsToImg(...
                size(alData.img),loc,params.spotSize);
        else
            spotsImg = zeros(size(alData.img));
            frameSize = size(alData.img);
            frameSize = frameSize(1:2);
            frameColNum = ismember(locVars,'frame');
            for i=1:alData.nFrames
                tmpLoc = loc(loc(:,frameColNum)==i,:);
                spotsImg(:,:,i) = addGaussianSpotsToImg(...
                    frameSize,tmpLoc,params.spotSize);
            end
        end
    end
end