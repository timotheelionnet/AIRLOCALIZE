function [loc,locVars] = run_gaussian_fit_on_all_spots_in_image3(...
    spotCandidates,alData,params)

%  loc contains the coordinates output by the fine localization algorithm (in pix)

tic
verbose = 0;

%% initialize arrays and set cutwidth: the 
nSpots = size(spotCandidates,1);       
switch params.fitMethod
    case '3DMaskFull'
        loc = zeros(nSpots,6+alData.isMovie);   
        locVars = {'x_in_pix','y_in_pix','z_in_pix',...
                'integratedIntensity','residuals','image_number'};
        if alData.isMovie
                locVars = [locVars,'frame'];
        end    
        cutWidth = [ceil(params.fittedRegionSize*params.psfSigma(1)),...
            ceil(params.fittedRegionSize*params.psfSigma(2))];        
    case '2DMaskOnLocalMaxProj'
        loc = zeros(nSpots,6+alData.isMovie); 
        locVars = {'x_in_pix','y_in_pix','z_in_pix',...
                'integratedIntensity','residuals','image_number'};
        if alData.isMovie
                locVars = [locVars,'frame'];
        end
        cutWidth = [ceil(params.fittedRegionSize*params.psfSigma(1)),...
            ceil(params.fittedRegionSize*params.psfSigma(2))];        
    case '3DGaussianFit'
        loc = zeros(nSpots,8+alData.isMovie);  
        locVars = {'x_in_pix','y_in_pix','z_in_pix',...
                'integratedIntensity','residuals','image_number'};
        if alData.isMovie
                locVars = [locVars,'frame'];
        end
        cutWidth = [ceil(params.fittedRegionSize*params.psfSigma(1)),...
            ceil(params.fittedRegionSize*params.psfSigma(2))];  
    case '2DGaussianMask'
        loc = zeros(nSpots,4+alData.isMovie);   
        locVars = {'x_in_pix','y_in_pix',...
                'integratedIntensity','residuals','image_number'};
        if alData.isMovie
                locVars = [locVars,'frame'];
        end
        cutWidth = ceil(params.fittedRegionSize*params.psfSigma(1));          
    case '2DGaussianFit'
        loc = zeros(nSpots,6+alData.isMovie);  
        locVars = {'x_in_pix','y_in_pix',...
                'integratedIntensity','residuals','image_number'};
        if alData.isMovie
                locVars = [locVars,'frame'];
        end
        cutWidth = ceil(params.fittedRegionSize*params.psfSigma(1)); 
    otherwise
        if ischar(params.fitMethod)
            disp(['fit entry : ''',params.fitMethod,...
                ''' is not a recognized string. ',...
                'It should be one of the following ',...
                '(making sure it matches the data dimensions):']);
        else
            disp(['fit entry : is not recognized. ',...
            'It should be one of the following ',...
            '(making sure it matches the data dimensions):']);
        end
        disp('    3DMaskFull');
        disp('    2DMaskOnLocalMaxProj');
        disp('    3DGaussianFit');
        disp('    2DGaussianMask');
        disp('    2DGaussianFit');     
end
     
%% loop over each pre-detected spot    

% pre-fill image and frame number
loc(:,params.numDim+3) = alData.fileIdx;  
if alData.isMovie
    loc(:,params.numDim+4) = alData.curFrame;  
end

for j=1:nSpots
    %progress update message
    if(mod(j,100)==0 && verbose==1) 
        disp(['  spot ',num2str(j),' out of ' num2str(nSpots)]); 
    end

    %Background correction
    switch params.bgCorrectionMode
        case 'localPlane' 
            [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                gen_linear_interpol_clean_small2(...
                alData,spotCandidates(j,1:params.numDim),cutWidth,...
                params.bgRegionThickness ,'large');
        case 'localMedian'
            [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                subtract_median2(...
                alData,spotCandidates(j,1:params.numDim),cutWidth,...
                params.bgRegionThickness ,'large','local');
        case 'globalMedian'
            [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                subtract_median2(...
                alData,spotCandidates(j,1:p.numdim),cutWidth,...
                params.bgRegionThickness ,'large','global');
    end
    
    %Gaussian Fit
    switch params.fitMethod
        case '3DMaskFull'
            % fine localization using the gaussian mask algorithm
            [x0,y0,z0,N0,err0] = gaussian_mask_small2(...
                stack_bg_corr,new_ctr,params);
            
            %storing the results
            loc(j,1) = x0 + ROIlimits(1,1) - 1;
            loc(j,2) = y0 + ROIlimits(1,2) - 1;
            loc(j,3) = z0 + ROIlimits(1,3) - 1;
            loc(j,4) = N0; %Itot
            loc(j,5) = err0; %Residuals        
        case '2DMaskOnLocalMaxProj'
            %local maxproj
            img_bg_corr = local_maxproj(...
                stack_bg_corr,new_ctr,params);

            % fine localization using the gaussian mask algorithm
            [x0,y0,~,N0,err0] = gaussian_mask_small2(...
                img_bg_corr,new_ctr(1:2),params);
            
            %storing the results
            loc(j,1) = x0 + ROIlimits(1,1) - 1;
            loc(j,2) = y0 + ROIlimits(1,2) - 1;
            loc(j,3) = new_ctr(3) + ROIlimits(1,3) - 1;
            loc(j,4) = N0;  %Itot
            loc(j,5) = err0;%Residuals
        case '3DGaussianFit'
            Gaussout = gaussian_fit_local3(stack_bg_corr,new_ctr,params,1);
            loc(j,1)=Gaussout(5) + ROIlimits(1,1) - 1;  %xc
            loc(j,2)=Gaussout(6) + ROIlimits(1,2) - 1;  %yc
            loc(j,3)=Gaussout(7) + ROIlimits(1,3) - 1;  %zc
            loc(j,4)=Gaussout(8); %Itot
            loc(j,5)=Gaussout(9); %Residuals
        case '2DGaussianMask'
            % fine localization using the gaussian mask algorithm
            [x0,y0,~,N0,err0] = gaussian_mask_small2(...
                stack_bg_corr,new_ctr,params);
            %storing the results
            loc(j,1)=x0 + ROIlimits(1,1) - 1;
            loc(j,2)=y0 + ROIlimits(1,2) - 1;
            loc(j,3)=N0;%Itot
            loc(j,4)=err0;%Residuals
        case '2DGaussianFit'
            Gaussout = gaussian_fit_local3(stack_bg_corr,new_ctr,params,1);
            loc(j,1)=Gaussout(4) + ROIlimits(1,1) - 1;  %xc
            loc(j,2)=Gaussout(5) + ROIlimits(1,2) - 1;  %yc
            loc(j,3)=Gaussout(6);  %Itot
            loc(j,4)=Gaussout(7);  %residuals       
    end
end 

%% array cleanup & report display
%setting the arrays to the right size, ordering them by decreasing intensity and removing the potential double identifications
[loc,nWrong,nDouble] = ...
    clean_up_spots_array_clean4(loc,size(alData.img),params); 

nspots_found = size(loc,1);

t = toc;  

disp(['  eliminated ',num2str(nWrong),' wrong / ',num2str(nDouble),...
    ' double identifications; there remain ',num2str(nspots_found),...
    ' spots on image ',num2str(alData.fileIdx)]);
disp(['  gaussian mask complete after ', num2str(t), ' s']);


end

function img = local_maxproj(alData,spot_ctr,params)

zc = spot_ctr(3);
cutwidth_z = ceil(params.psfSigma(2) * params.fittedRegionSize);
[~,~,nz]= size(alData.img);
    
%the number of planes max projected
ROIlimits(1,3) = ceil(zc) - floor(cutwidth_z);
ROIlimits(1,3) = max(1,ROIlimits(1,3));
zmax = ceil(zc) + floor(cutwidth_z);
zmax = min(nz,zmax);

img = max( alData.img(:,:,ROIlimits(1,3):zmax),[],3 );
end

