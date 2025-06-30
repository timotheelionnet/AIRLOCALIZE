function params = explore_stack_2channels4(params,alData)

%% creating figure
fh = figure(...
              'Units','characters',...
              'MenuBar','none',...
              'Toolbar','none',...
              'NumberTitle','off',...
              'Position',[20 3 200 55],...
              'Visible','off'); 
          
%% parsing arguments 
%making sure that raw data and smoothed image are loaded into the
%airlocalizeData object
overwrite = 1;
if isempty(alData.img)
    alData.retrieveImg(overwrite);
    alData.retrieveSmoothedImg(overwrite);
elseif isempty(alData.smooth)
    alData.retrieveSmoothedImg(overwrite);
end

if ~ismember(ndims(alData.img),[2,3]) ...
        || ~ismember(ndims(alData.smooth),[2,3]) 
    disp(['input data size is', num2str(size(alData.img)),...
        ' it should be 2D or 3D.']);
    return;
elseif ndims(alData.img) ~= ndims(alData.smooth)
    disp(['input data size is', num2str(size(alData.img)),...
        ' it should be 2D or 3D.']);
end

[nx, ny, nz] = size(alData.img);
setappdata(fh,'nx',nx);
setappdata(fh,'ny',ny);
setappdata(fh,'nz',nz);
xfreeze = round(nx/2); 
yfreeze = round(ny/2); 
setappdata(fh,'xfreeze',xfreeze);
setappdata(fh,'yfreeze',yfreeze);
setappdata(fh,'hIminVal',min(alData.img(:)));
setappdata(fh,'hImaxVal',max(alData.img(:)));
setappdata(fh,'hthresh_levelVal',params.threshLevel);
thresh = struct();
thresh.level = params.threshLevel;
thresh.units = params.threshUnits;
setappdata(fh,'defaultBackgroundColor','k');

if strcmp(thresh.units,'legacySD')
    thresh.sd = std(alData.smooth(:,:,1),0,'all'); % using first frame to compute sd
else
    thresh.sd = std(alData.img(:,:,1),0,'all'); % using first frame to compute sd
end
setappdata(fh,'thresh',thresh);

fh = set_figure_and_panel_for_explore_stack_set_threshold(...
    alData,fh,params);
handles = guihandles(fh);


%%%%%%%%%%%%%%%%%%%%%%%%%%% set all the callback functions %%%%%%%%%%%%%%%%%            
set(fh,'WindowButtonMotionFcn', {@viewlocaldata,alData,fh});
set(fh,'WindowButtonUpFcn', {@change_state,fh});   
set(handles.hclose,'Callback', {@clean_close,fh});   
set(handles.zslider,'Callback',{@slider_change,alData,fh});
set(handles.hzpos,'Callback',{@zpos_change,alData,fh});            
set(handles.hImin,'Callback',{@change_Iminmax,alData,fh});
set(handles.hImax,'Callback',{@change_Iminmax,alData,fh});
set(handles.hslice_auto,'Callback',{@slice_auto_adjust,alData,fh});
set(handles.hglob_auto,'Callback',{@global_auto_adjust,alData,fh});
set(handles.hxplot_range,'Callback',{@change_local_range,alData,fh});
set(handles.hyplot_range,'Callback',{@change_local_range,alData,fh});
set(handles.hzplot_range,'Callback',{@change_local_range,alData,fh});
set(handles.hoverlay,'Callback',{@toggle_overlay,fh,alData});
set(handles.hthresh_level,'Callback',{@update_threshold,fh,alData,params});
set(handles.hthresh_units,'SelectionChangeFcn',{@update_threshold,fh,alData,params});
set(fh,'WindowKeyPressFcn',{@moveCursorWithArrows,alData,fh});

set(fh,'SelectionType','alt');
set(fh,'CurrentAxes',handles.ha);
    
plot_current_z_stack_two_channels(fh,alData,1);
set(fh,'Visible','on');

uiwait;

thresh = getappdata(fh,'thresh');
params.threshLevel = thresh.level;
params.threshUnits = thresh.units;

close(fh);
drawnow;

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function moveCursorWithArrows(src,eventData,alData,fh)
    % list of keys that when pressed, map to actions
    keyList = {'uparrow','downarrow','leftarrow','rightarrow','h'};
    charList = {'<', ',','>', '.'};
    charKeyMap = {'comma','comma','period','period'};
    if ~ismember(eventData.Key,keyList) ...
            && ~ismember(eventData.Character,charList)
        return
    else
        if sum(ismember(keyList,eventData.Key)) == 1
            curKey = keyList{ismember(keyList,eventData.Key)};
        else
            curKey = charKeyMap{ismember(charList,eventData.Character)};
        end
    end

    handles = guihandles(fh);
    % key 'h' is the shortcut to toggle the overlay
    if strcmp(curKey,'h')
        set(handles.hoverlay,'Value', 1 - get(handles.hoverlay,'Value'));
        toggle_overlay([],[],fh,alData);
        return
    end
    
    if ndims(alData.img)==3
        [nx, ny, nz] = size(alData.img);
    else
        [nx, ny] = size(alData.img);
    end

    val = get(handles.hstate,'Value');
    if val == 2
        % in snap mode, move the z,y with arrows for fine positioning of
        % the cursor (useful for centering on the spots for gaussian fits)
        % move z plane with < and >
        xfreeze = getappdata(fh,'xfreeze');
        yfreeze = getappdata(fh,'yfreeze');
        z = str2double(get(handles.hzpos,'String'));
        z = round(z);
        switch curKey
            case 'uparrow'
                xfreeze = xfreeze-1;
            case 'downarrow'
                xfreeze = xfreeze+1;
            case 'leftarrow'
                yfreeze = yfreeze-1;
            case 'rightarrow'
                yfreeze = yfreeze+1;
            case {'comma'}
                z = z - 1;
            case {'period'}  
                z = z + 1;
        end

        % enforce that the new cursor position remains within bounds
        xfreeze = max(min(xfreeze,nx),1);
        yfreeze = max(min(yfreeze,ny),1);
        if ndims(alData.img)==3
            z = max(min(z,nz),1);
            set(handles.zslider,'Value',z);
            set(handles.hzpos,'String',num2str(z));
            slider_change([],[],alData,fh);
        else 
            z = 1;
        end
        
        % update the figure data with the new cursor position
        setappdata(fh,'xfreeze',xfreeze);
        setappdata(fh,'yfreeze',yfreeze);
        set(handles.hxpos,'String',num2str(xfreeze));
        set(handles.hypos,'String',num2str(yfreeze));
        Int = alData.img(xfreeze,yfreeze,z);
        set(handles.hIval,'String',num2str(Int));

        plot_local_data_two_channels(xfreeze,yfreeze,z,fh,alData);
    else
        % in grab mode, move the z-slider if the image is 3D
        if ndims(alData.img)==3
            z = str2double(get(handles.hzpos,'String'));
            z = round(z);
            switch curKey
                case {'leftarrow','comma'}
                    if z > 1
                        z = z-1;
                        set(handles.zslider,'Value',z);
                        set(handles.hzpos,'String',num2str(z));
                        slider_change([],[],alData,fh);
                    end
                case {'rightarrow','period'}
                    if z<nz
                        z = z+1;
                        set(handles.zslider,'Value',z);
                        set(handles.hzpos,'String',num2str(z));
                        slider_change([],[],alData,fh);
                    end
            end
        end
    end
end

function update_threshold(src,eventdata,fh,alData,params)
    handles = guihandles(fh);
    
    % update the threshold value from GUI
    thresh = getappdata(fh,'thresh');
    thresh.level = getCleanNumberFromEditEntry(handles,'hthresh_level',fh,'hthresh_levelVal');
    % update the threshold units from GUI
    str_units = get(get(handles.hthresh_units,'SelectedObject'),'Tag');
    % if user selects adaptive but this wasnt set up in the parameters,
    %this will lead to problems bc there is no mask to use; revert to SD
    if strcmp(str_units,'adaptive') && ~strcmp(params.threshUnits,'adaptive')
        str_units = 'SD';
        allButtons = handles.hthresh_units.Children;
        sdButton = findobj(allButtons,'Tag','SD');        
        set(handles.hthresh_units,'SelectedObject',sdButton);
        %drawnow;
    end
    
    switchThresholdImage = 0;
    if strcmp(params.threshUnits,'adaptive')        
        if strcmp(str_units,'adaptive') ...
                && ismember(thresh.units,{'absolute','SD','legacySD'}) ...
                || strcmp(thresh.units,'adaptive') ...
                && ismember(str_units,{'absolute','SD','legacySD'}) 
            switchThresholdImage = 1; 
        end
    end
    if strcmp(str_units,'absolute') 
        thresh.units = 'absolute';
    elseif strcmp(str_units,'SD') 
        thresh.units = 'SD';
    elseif strcmp(str_units,'adaptive') 
        thresh.units = 'adaptive';
    elseif strcmp(str_units,'legacySD') 
        thresh.units = 'legacySD';
    end
    % if we switch between adaptive and classic thresholding,
    % we need to re-generate the smoothed image.
    if switchThresholdImage 
        set(handles.hthresh_units,'Title','Loading Smoothed Image...');
        drawnow;
        params2 = params.duplicate();
        params2.threshUnits = thresh.units;
        overwrite = 1;
        alData.retrieveSmoothedImg(params2,overwrite);
        set(handles.hthresh_units,'Title','Units');
        drawnow;
    end
    
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    
    % update current SD value that will be used to plot selected pixels
    if strcmp(thresh.units,'SD')
        if ndims(alData.img) == 3 && ~alData.isMovie
            thresh.sd = std(alData.img,0,'all');
        else
            thresh.sd = std(alData.img(:,:,z),0,'all');
        end
    elseif strcmp(thresh.units,'legacySD')
        if ndims(alData.img) == 3 && ~alData.isMovie
            thresh.sd = std(alData.smooth,0,'all');
        else
            thresh.sd = std(alData.smooth(:,:,z),0,'all');
        end
    end
    
    % save new thresh settings to variable
    setappdata(fh,'thresh',thresh);
    
    % compute number of voxels above threshold
    nVox = computeNumberOfVoxelsAboveThreshold(thresh,alData);
    
    set( handles.hnpix_above_threshold,'String', num2str(nVox));
    
    %replot the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_local_data_two_channels(x,y,z,fh,alData);
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);
    
    %cleaning the global link 
    clear('thresh','stack2','xfreeze','yfreeze','handles','thresh',...
        'x','y','z','npix_txt','npix_above_threshold');
end

function nVox = computeNumberOfVoxelsAboveThreshold(thresh,alData)
    if strcmp(thresh.units,'SD')
        if ndims(alData.img) == 3 && ~alData.isMovie
            threshInt = thresh.level * std(alData.img,0,'all');
        else
            threshInt = thresh.level * ones(size(alData.smooth)) ...
                .* repmat(std(alData.img,0,[1,2]),...
                size(alData.smooth,1),size(alData.smooth,2),1);
        end
    elseif strcmp(thresh.units,'legacySD')
        if ndims(alData.img) == 3 && ~alData.isMovie
            threshInt = thresh.level * std(alData.smooth,0,'all');
        else
            threshInt = thresh.level * ones(size(alData.smooth)) ...
                .* repmat(std(alData.smooth,0,[1,2]),...
                size(alData.smooth,1),size(alData.smooth,2),1);
        end
    elseif strcmp(thresh.units,'absolute') || strcmp(thresh.units,'adaptive')
        % note that in adaptive mode the smoothed image is pre-normalized
        % within individual ROIs so we are just applying the same threshold
        % across the board
        threshInt = thresh.level;
    end
    
    nVox = sum(alData.smooth > threshInt,'all');

end

function toggle_overlay(src,eventdata,fh,alData)
    handles = guihandles(fh);   
    ny = size(alData.img,2);
    
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    
    %switching the state of the button
    is_overlay_on = 1 - get(handles.hoverlay,'Value');
    if is_overlay_on == 1 
        set(handles.hoverlay,'String','Hide Overlay (shortcut:h)');
    else 
        set(handles.hoverlay,'String','Show Overlay (shortcut:h)');
    end
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    
    val = get(handles.hstate,'Value');  %this is the snap/grab mode for the local data panel
    if val == 1                 %if grab mode is activated, I replot the zoom data at some arbitrary location
        %plot_zoom_two_channels(1,ny,z,zoomwidth,alData,fh);
        viewlocaldata([],[],alData,fh);
          
    else                %if snap mode is activated, I plot the local data at the same location
        plot_zoom_two_channels(...
            getappdata(fh,'xfreeze'),getappdata(fh,'yfreeze'),z,...
        zoomwidth,alData,fh);
    end
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);
    
    clear('xfreeze','yfreeze','stack2','handles',...
        'nx','ny','nz','is_overlay_on','z','val','zoomwidth');
end

function plot_current_z_stack_two_channels(fh,alData,z)

    thresh = getappdata(fh,'thresh');
    if strcmp(thresh.units,'absolute') || strcmp(thresh.units,'adaptive')
        threshInt = thresh.level;
    elseif strcmp(thresh.units,'SD') || strcmp(thresh.units,'legacySD')
        threshInt = thresh.level*thresh.sd;
    end
    
    handles = guihandles(fh);    
    %% parsing dimensions of input stack   
    nd = ndims(alData.img);
    if nd == 3
        [nx, ny, nz] = size(alData.img);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
    elseif nd == 2
        [nx, ny] = size(alData.img);
    else
        disp('wrong stack size');
    end
      
    Imin = getCleanNumberFromEditEntry(handles,'hImin',fh,'hIminVal');
    Imax = getCleanNumberFromEditEntry(handles,'hImax',fh,'hImaxVal');
    
    %% filling data 
    curSlice = zeros(nx,ny,3);
    
    %red channel: overlay or nothing
    is_overlay_on = 1 - get(handles.hoverlay,'Value');
    if ~is_overlay_on
        threshInt = Inf;
    end
    
    if nd == 3 
        curSlice(:,:,1) = alData.smooth(:,:,z)>threshInt;
    else
        curSlice(:,:,1) = alData.smooth(:,:)>threshInt;
    end
    
    %green channel: data
    if(Imin >= Imax)
        curSlice(:,:,2) = 0;
    else
        if nd == 3
            curSlice(:,:,2) = ...
                min( max(( alData.img(:,:,z) -Imin)/(Imax-Imin),0) , 1);
        else
            curSlice(:,:,2) = ...
                min( max(( alData.img(:,:) -Imin)/(Imax-Imin),0) , 1);
        end
    end
    
    %nothing in the blue channel
    curSlice(:,:,3) = 0; 
    
    %%    
    ha = handles.ha;
    imagesc(ha,curSlice); 
    xlim(ha,[1,max(nx,ny)]);
    ylim(ha,[1,max(nx,ny)]);
    defaultBackgroundColor = getappdata(fh,'defaultBackgroundColor');
    set(ha,'Color',defaultBackgroundColor);
    set(ha,'Tag','ha');
        
    clear('cur_slice','nx','ny','nz','Imin','Imax','newXLim','newYLim','stack1');
    clear('stack2','handles','ha','thresh','threshInt','is_overlay_on');
end

function plot_zoom_two_channels(x,y,z,width,alData,fh)

 %inputs a 3D or 2D image and plots a locally  zoomed version  
    handles = guihandles(fh);
    zoomh = handles.zoomh;
    
    thresh = getappdata(fh,'thresh');
    if strcmp(thresh.units,'absolute') || strcmp(thresh.units,'adaptive')
        threshInt = thresh.level;
    elseif strcmp(thresh.units,'SD') || strcmp(thresh.units,'legacySD')
        threshInt = thresh.level*thresh.sd;
    end
    nd = ndims(alData.img);
    if nd == 3
        [nx, ny, nz] = size(alData.img);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end   
    elseif nd ==2
        [nx, ny] = size(alData.img);
    end
    
    Imin = getCleanNumberFromEditEntry(handles,'hImin',fh,'hIminVal');
    Imax = getCleanNumberFromEditEntry(handles,'hImax',fh,'hImaxVal');
    
%% computing the value of the local slice 
    
    curSlice = zeros(nx,ny,3);

    %red channel: overlay
    is_overlay_on = 1 - get(handles.hoverlay,'Value');
    if ~is_overlay_on 
        threshInt = Inf;
    end
    
    if nd == 2
        curSlice(:,:,1) = alData.smooth > threshInt;
    else
        curSlice(:,:,1) = alData.smooth(:,:,z) > threshInt;
    end

    %green channel: data
    if(Imin >= Imax)
        curSlice(:,:,2) = 0;
    else
        if nd == 2
            curSlice(:,:,2) = ...
                min(max( ( alData.img-Imin) / (Imax-Imin),0),1);
        else
            curSlice(:,:,2) = ...
                min(max( ( alData.img(:,:,z) -Imin) / (Imax-Imin),0),1);
        end
    end

    %nothing in the blue channel
    curSlice(:,:,3) = 0;      
    
%% chosing the window corresponding to the cursor position
    if 2*width+1>nx
        if 2*width+1>ny
            image(zoomh,curSlice);
        elseif (y - width >= 1) && (y + width <= ny)
            image(zoomh,[y-width y+width],[1 nx],...
                curSlice(1:nx,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image(zoomh,[yy-width,yy+width],[1,nx],...
                curSlice(1:nx,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image(zoomh,[yy-width,yy+width],[1,nx],...
                curSlice(1:nx,yy-width:yy+width,1:3));
        end
    elseif (x - width >= 1) && (x + width <= nx)
        if 2*width+1>ny
            image(zoomh,[1,ny],[x-width,x+width],...
                curSlice(x-width:x+width,1:ny,1:3));
        elseif (y - width >= 1) && (y + width <= ny)
            image(zoomh,[y-width,y+width],[x-width,x+width],...
                curSlice(x-width:x+width,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image(zoomh,[yy-width,yy+width],[x-width,x+width],...
                curSlice(x-width:x+width,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image(zoomh,[yy-width,yy+width],[x-width,x+width],...
                curSlice(x-width:x+width,yy-width:yy+width,1:3));
        end
    elseif (x - width < 1) && (x + width <= nx) 
        xx = width+1;
        if 2*width+1>ny
            image(zoomh,[1,ny],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,1:ny,1:3));
        elseif (y - width >= 1) && (y + width <= ny)
            image(zoomh,[y-width,y+width],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image(zoomh,[yy-width,yy+width],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image(zoomh,[yy-width,yy+width],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,yy-width:yy+width,1:3));
        end
    elseif (x + width > nx) && (x - width >= 1)
        xx = nx - width;
        if 2*width+1>ny
            image(zoomh,[1,ny],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,1:ny,1:3));
        elseif (y - width >= 1) && (y + width <= ny)
            image(zoomh,[y-width,y+width],[xx-width,xx+width],...
            curSlice(xx-width:xx+width,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image(zoomh,[yy-width,yy+width],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image(zoomh,[yy-width,yy+width],[xx-width,xx+width],...
                curSlice(xx-width:xx+width,yy-width:yy+width,1:3));
        end
    end 
    set(zoomh,'Tag','zoomh');
    clear('cur_slice','nx','ny','nz','Imin','Imax','is_overlay_on','stack1','stack2');
    clear('xx','yy','handles','thresh','threshInt','zoomh');
    clear('x','y','z','width','zoomh','stack1','stack2','fh','hImin','hImax','hoverlay');
end

function nOut = getCleanNumberFromEditEntry(handles,tagName,fh,valName)
    % reads the String field of a 'text' object with 'tagName' from the
        % handles 'handles.
        % cleans up the string and converts to a numeral nOut. Sets the
        % appdata property 'valName' of fh to nOut.
        % if nOut is NaN, reverts nOut to the value stored in the 'valName'
        % property of fh.

    strIn = get(handles.(tagName),'String');
    
    % ensures any non-numeric characters that were entered accidentally via hotkeys 
    % are eliminated
    strIn = regexprep(strIn, '[^0-9.]', '');
    
    % remove all dots except the leftmost one 
    % Find indices of all dots
    dotIdx = strfind(strIn, '.');
    
    % If more than one dot, remove all except the first
    if numel(dotIdx) > 1
        remove = dotIdx(2:end);
    elseif isscalar(dotIdx)
        if dotIdx(1) == length(strIn)
            remove = dotIdx;
        else
            remove = [];
        end
    else
        remove = [];
    end
    strIn(remove) = [];  % remove unwanted dots

    set(handles.(tagName),'String',strIn); % clean up entry if needed
    strIn = get(handles.(tagName),'String');
    % Set selection to the end of the string
    nOut = str2double(strIn);
    if ~isnan(nOut)
        setappdata(fh,valName,nOut);
    else
        % revert to previously vetted value of the entry
        nOut = getappdata(fh,valName);
        set(handles.(tagName),'String',num2str(nOut));
    end
end

function change_local_range(src,eventdata,alData,fh)
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_local_data_two_channels(x,y,z,fh,alData);
    clear('x','y','z','handles');
end

function clean_close(src,eventdata,fh)
    uiresume;
end

function change_state(src,eventdata,fh)
    %switches between the snap and grab mode
    %hstate value = 1: grab
    %hstate value = 2: snap
    
    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    % check the current state; 2 = snap; 1 = grap
    val = get(handles.hstate,'Value');
    val = 3-val; %1->2 or 2->1
    set(handles.hstate,'Value',val);
    
    if val == 2   %if switching to snap state
        set(fh,'CurrentAxes',handles.ha);
        point = get(gca,'CurrentPoint');
        yfreeze = round(point(1,1));    %inverse convention for mouse and array (xy -> yx)
        xfreeze = round(point(1,2));
        xfreeze = max(1,min(xfreeze,nx));
        yfreeze = max(1,min(yfreeze,ny));
        
    end
    setappdata(fh,'xfreeze',xfreeze);
    setappdata(fh,'yfreeze',yfreeze);
    
    clear('point','val','point','xfreeze','yfreeze');
    clear('nx','ny','handles');
end

function slider_change(zslider,eventdata,alData,fh)
    nd = ndims(alData.img);
    if nd == 2
        return
    end
    
    handles = guihandles(fh);
    
    % updating the z-slider and the z info window
    z = get(handles.zslider,'Value');
    z = round(z);
    set(handles.zslider,'Value',z);
    set(handles.hzpos,'String',num2str(z));
    
    % updating the Imin/Imax info
    [curImin,curImax,curImed,curIsd] = getStats(alData,z);
    
    set(handles.hslice_Imin_Info_val,'String',num2str(curImin));
    set(handles.hslice_Imed_Info_val,'String',num2str(curImed));
    set(handles.hslice_Imax_Info_val,'String',num2str(curImax)); 
    set(handles.hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    % if movie, updating the current thresh.sd value
    if alData.isMovie
        thresh = getappdata(fh,'thresh');
        if strcmp(thresh.units,'legacySD')
            thresh.sd = std(alData.smooth(:,:,z),0,'all');
        elseif strcmp(thresh.units,'SD')
            thresh.sd = std(alData.img(:,:,z),0,'all');
        end
        setappdata(fh,'thresh',thresh);
    end
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    if nd == 3
        Int = alData.img(x,y,z);
    elseif nd == 2
        Int = alData.img(x,y);
    end
    set(handles.hIval,'String',num2str(Int));

    plot_local_data_two_channels(x,y,z,fh,alData);
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);  
    
    clear('x','y','z','Int','curImin','curImax','curImed','curIsd');
    clear('nx','ny','xfreeze','yfreeze','handles');
end

function [curImin,curImax,curImed,curIsd] = getStats(alData,z)
    nd = ndims(alData.img);
    if nd == 3
        curImin = min(alData.img(:,:,z),[],'all'); 
        curImax = max(alData.img(:,:,z),[],'all'); 
        curImed = median(alData.img(:,:,z),'all'); 
        curIsd = std(alData.img(:,:,z),0,'all');
    else
        curImin = min(alData.img(:,:),[],'all'); 
        curImax = max(alData.img(:,:),[],'all'); 
        curImed = median(alData.img(:,:),'all'); 
        curIsd = std(alData.img(:,:),0,'all');
    end
end

function zpos_change(hzpos,eventdata,alData,fh)
    nd = ndims(alData.img);
    if nd == 2
        return
    end
    nz = getappdata(fh,'nz');
    
    handles = guihandles(fh);
    
    %updating the z-slider and the z info window
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    z = max(min(z,nz) , 1);
    set(handles.hzpos,'String',num2str(z));
    set(handles.zslider,'Value',z);
    
    %updating the Imin/Imax info
    [curImin,curImax,curImed,curIsd] = getStats(alData,z);
    
    set(handles.hslice_Imin_Info_val,'String',num2str(curImin));
    set(handles.hslice_Imed_Info_val,'String',num2str(curImed));
    set(handles.hslice_Imax_Info_val,'String',num2str(curImax));
    set(handles.hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    % if movie, updating the current thresh.sd value
    if alData.isMovie
        thresh = getappdata(fh,'thresh');
        if strcmp(thresh.units,'legacySD')
            thresh.sd = std(alData.smooth(:,:,z),0,'all');
        elseif strcmp(thresh.units,'SD')
            thresh.sd = std(alData.img(:,:,z),0,'all');
        end
        setappdata(fh,'thresh',thresh);
    end
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    if nd == 3
        Int = alData.img(x,y,z);
    elseif nd == 2
        Int = alData.img(x,y);
    end
    set(handles.hIval,'String',num2str(Int));

    plot_local_data_two_channels(x,y,z,fh,alData);
    
    clear('x','y','z','Int','curImin','curImax','curImed','curIsd'); 
    clear('nx','ny','xfreeze','yfreeze','handles');
end

function change_Iminmax(hObj,eventdata,alData,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_zoom_two_channels(x,y,z,zoomwidth,alData,fh);
    
    clear('x','y','z','zoomwidth');
    clear('handles');
end

function slice_auto_adjust(hslice_auto,eventdata,alData,fh)
    nd = ndims(alData.img);
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %updating Imin / Imax info
    [Imin,Imax,~,~] = getStats(alData,z);
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_two_channels(x,y,z,zoomwidth,alData,fh);
    clear('x','y','z','zoomwidth','Imin','Imax');
    clear('handles');
end

function global_auto_adjust(hglob_auto,eventdata,alData,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %updating the Imin/Imax info
    Imin = min(alData.img,[],'all'); 
    Imax = max(alData.img,[],'all');
    
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,alData,z);
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_two_channels(x,y,z,zoomwidth,alData,fh);
    
    clear('x','y','z','zoomwidth','Imin','Imax');
    clear('handles');
end

function viewlocaldata(src,eventdata,alData,fh)
    %if checkIfCallBackIsTooFrequent(src,t,minT) || isMultipleCall
    if isMultipleCall
        return
    end
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    val = get(handles.hstate,'Value');          %this is the snap/grab mode
    ha = handles.ha;
    xh = handles.xh;
    yh = handles.yh;
    zh = handles.zh;
    zoomh = handles.zoomh;
    
    %I replot things only if the grab mode is active
    if val == 1 

        %I get the current (x,y,z) 
        [x,y] = get_current_pointer_position(fh); 
        z = get(handles.zslider,'Value');
        z = round(z);

        %I update the current display of the pointer info
        set(handles.hxpos,'String',num2str(x));
        set(handles.hypos,'String',num2str(y));

        if(x>=1 && y >=1 && x <= nx && y <= ny)
            if ndims(alData.img) == 3
                Int = alData.img(x,y,z);
            else
                 Int = alData.img(x,y);
            end
            set(handles.hIval,'String',num2str(Int));

            %plotting the profiles and zoom window
            plot_local_data_two_channels(x,y,z,fh,alData);
        end
    end
    
    % this command is key to ensure that there is no pile up of callbacks
    % which leads to lags
    drawnow;
    
    % reset the axes tags because plotting messes with them
    set(ha,'Tag','ha');
    set(xh,'Tag','xh');
    set(yh,'Tag','yh');
    set(zh,'Tag','zh');
    set(zoomh,'Tag','zoomh');
    
    clear('x','y','z','val','Int','nx','ny','handles');
end

function plot_local_data_two_channels(x,y,z,fh,alData)
   
    handles = guihandles(fh);
    
    xrange = str2double(get(handles.hxplot_range,'String'));
    yrange = str2double(get(handles.hyplot_range,'String'));
    if ndims(alData.img)==3
        zrange = str2double(get(handles.hzplot_range,'String'));
    end
    zoomwidth = str2double(get(handles.hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile(x,y,z,xrange,alData,fh);
    plot_yprofile(x,y,z,yrange,alData,fh);
    if ndims(alData.img)==3
        plot_zprofile(x,y,z,zrange,alData,fh);
    end
    plot_zoom_two_channels(x,y,z,zoomwidth,alData,fh);
    
    clear('handles','xrange','yrange','zrange','zoomwidth');
end

function [xarr, yarr] = get_current_pointer_position(fh)

    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    
    val = get(handles.hstate,'Value');%this is the snap/grab mode for the local data panel

    if val == 1         %if grab mode is activated, I replot the local data at some arbitrary location
            set(fh,'CurrentAxes',handles.ha);
            point = get(gca,'CurrentPoint');
            xpic = round(point(1,1));   
            ypic = round(point(1,2));
            xarr = ypic;    
            yarr = xpic;    %pic and array have different xy conventions
            xarr = max(xarr,1);  
            xarr = min(xarr,nx);
            yarr = max(yarr,1);  
            yarr = min(yarr,ny);
    else
            xarr = xfreeze; 
            yarr = yfreeze;
    end

    clear('xfreeze','yfreeze','nx','ny','handles','xpic','ypic','val');
end

function plot_xprofile(xarr,yarr,z,xrange,alData,fh)
    
    handles = guihandles(fh);
    xh = handles.xh;
    
    nx = size(alData.img,1);
    if ndims(alData.img) == 3
        xprofile(1:nx) = alData.img(1:nx,yarr,z);
    else
        xprofile(1:nx) = alData.img(1:nx,yarr);
    end

    x1 = max(1,xarr-xrange);
    x2 = min(nx,xarr+xrange);
    Imin = min(xprofile(x1:x2));
    Imax = max(xprofile(x1:x2));
    
    plot(xh,x1:x2,xprofile(x1:x2),[xarr, xarr],[Imin, Imax]);
    xlim(xh,[x1,x2]);
    set(xh,'Tag','xh');
    clear('nx','xprofile','Imin','Imax','x','handles','xh');
end

function plot_yprofile(xarr,yarr,z,yrange,alData,fh)

    handles = guihandles(fh);
    yh = handles.yh;
    ny = size(alData.img,2);
    if ndims(alData.img) == 3
        yprofile(1:ny) = alData.img(xarr,1:ny,z);
    else
        yprofile(1:ny) = alData.img(xarr,1:ny);
    end
    
    y1 = max(1,yarr-yrange);
    y2 = min(ny,yarr+yrange);
    Imin = min(yprofile(y1:y2));
    Imax = max(yprofile(y1:y2));
    
    plot(yh,y1:y2,yprofile(y1:y2),[yarr, yarr],[Imin, Imax]); 
    xlim(yh,[y1,y2]);
    set(yh,'Tag','yh');
    clear('ny','yprofile','Imin','Imax','y','handles','yh');
end

function plot_zprofile(xarr,yarr,z,zrange,alData,fh)
    if ndims(alData.img) ~=3
        return
    end
    handles = guihandles(fh);
    zh = handles.zh;
    
    nz = size(alData.img,3);
    zprofile(1:nz) = alData.img(xarr,yarr,1:nz);
    
    z1 = max(1,z-zrange);
    z2 = min(nz,z+zrange);
    Imin = min(zprofile(z1:z2));
    Imax = max(zprofile(z1:z2));
    
    plot(zh, z1:z2, zprofile(z1:z2),[z, z],[Imin, Imax]); 
    xlim(zh,[z1,z2]);
    set(zh,'Tag','zh');
    clear('nz','zprofile','Imin','Imax','z','zz','handles','zh');
end

function flag = isMultipleCall()
    s = dbstack();
    % s(1) corresponds to isMultipleCall
    if numel(s)<=2 
        flag=false; 
        return; 
    end
    % compare all functions on stack to name of caller
    count = sum(strcmp(s(2).name,{s(:).name}));
    % is caller re-entrant?
    if count>1
        flag=true; 
    else
        flag=false; 
    end
end

function fh = set_figure_and_panel_for_explore_stack_set_threshold(...
    alData,fh,params)


%%%%%%%%%%%%%%%%%%%%% variables initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nd = ndims(alData.img);
if nd == 3
    [~,~,nz] = size(alData.img);
else
    nz = 1;
end

strXpos = '1';     % default x value of the current pixel
strYpos = '1';     % default y value of the current pixel
strIval = num2str(alData.img(1,1,1));   % value of the intensity @ the current pixel
strZpos = '1';     % the default value for the z stack displayed

%values of global min max and median
Imed = median(alData.img,'all');
Isd = std(alData.img,0,'all');
Imin = min(alData.img,[],'all');
Imax = max(alData.img,[],'all');

strImin = num2str(Imin,'%10.4f');
strImax = num2str(Imax,'%10.4f');
strImed = num2str(Imed,'%10.4f');
strIsd = num2str(Isd,'%10.4f');

% values of min max and median of current slice 
z = 1; % (z=1 when you open the viewer)
[curImin,curImax,curImed,curIsd] = getStats(alData,z);

strcurImin = num2str(curImin,'%10.4f');
strcurImax = num2str(curImax,'%10.4f');
strcurImed = num2str(curImed,'%10.4f');
strcurIsd = num2str(curIsd,'%10.4f');

%default values of the ranges of the 3 axial views
xpix = num2str(50);
ypix = num2str(50);
zpix = num2str(50);

%default value of the size of the zoom window.
startzoomwidth = num2str(50);

%% %%%%%%%%%%%%%%%%%%%%%%% main figure, plots and title %%%%%%%%%%%%%%%%%%%%%%%  

%set(fh,'Toolbar','figure');
set(fh,'Name',alData.curImgFile);              
axes('Parent',fh,'Units','characters','Color', 'k',...
    'Position',[10,6.15,102.4,39.4],'Tag','ha');     %main figure
        
axes('Parent',fh,'Units','characters',...
    'Position',[120 45.4 66 7.7],'Tag','xh'); 
axes('Parent',fh,'Units','characters',...
    'Position',[120 33.8 66 7.7],'Tag','yh'); 
az = axes('Parent',fh,'Units','characters',...
    'Position',[120 22.3 66 7.7],'Tag','zh'); 
if nd == 2
    set(az,'visible','off');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% settings for the axial views %%%%%%%%%%%%%%%%
uicontrol('Units','characters',...
                'Style','text','String','profile along x',...
                'Position',[144,53.3,20,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','profile along y',...
                'Position',[144,41.8,20,1.2]);
            
uprofstrZ = uicontrol('Units','characters',...
                'Style','text','String','profile along z',...
                'Position',[144,30.3,20,1.2]);
if nd == 2
    set(uprofstrZ,'visible','off');
end

uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,49.2,10,3.6]); 
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,37.7,10,3.6]); 
unpixZstr = uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,26.2,10,3.6]); 
if nd == 2
    set(unpixZstr,'visible','off');
end

uicontrol('Units','characters',...
                'Tag','hxplot_range',...
                'Style','edit','String',xpix,...
                'Position',[190,47.7,6,1.2]);                  
uicontrol('Units','characters',...
                'Tag','hyplot_range',...
                'Style','edit','String',ypix,...
                'Position',[190,36.2,6,1.2]);
            
uzpix = uicontrol('Units','characters',...
                'Tag','hzplot_range',...
                'Style','edit','String',zpix,...
                'Position',[190,24.6,6,1.2]);
if nd == 2
    set(uzpix,'visible','off');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% zoom plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
axes('Parent',fh,'Units','characters',...
    'Position',[144 1.9 46 17.7],'Tag','zoomh');             

uicontrol('Units','characters',...
                'Style','text','String','window size',...
                'Position',[192,12,8,2.3]);
            
uicontrol('Units','characters',...
                'Tag','hzoom_width',...
                'Style','edit','String',startzoomwidth,...
                'Position',[192,10.4,6,1.2]); 

%% %%%%%%%%%%%%%%%%%%%%%%%% Threshold Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hthresh = uipanel('parent',fh,'Title','Detection Threshold',...
    'Units','characters',...
    'Position',[116,4,25,12]);

uicontrol('parent',hthresh,'Units','characters',...
            'Style','text','String','Value',...
            'Position',[1,9,8,1.2]);   
        
uicontrol('parent',hthresh,'Units','characters',...
            'Tag','hthresh_level',...
            'Style','edit','String',num2str(params.threshLevel),...
            'Position',[11,9,8,1.2]);
                
thresh = getappdata(fh,'thresh');   
nVox = computeNumberOfVoxelsAboveThreshold(thresh,alData);
                    
uicontrol('parent',hthresh,'Units','characters',...
            'Style','text','String','# of Voxels',...
            'Position',[1,7,15,1.2]);     
        
uicontrol('parent',hthresh,'Units','characters',...
            'Tag','hnpix_above_threshold',...
            'Style','text','String',num2str(nVox),...
            'Position',[16,7,8,1.2]); 
        

hth = uibuttongroup('parent',hthresh,'Units','characters',...
    'Tag','hthresh_units',...
    'Enable','on',...
    'Title','Units',...
    'Position',[0.5 0.2 20 7]);

hAbs = uicontrol('parent',hth,'Style','Radio','Units','characters',...
        'String','Absolute',...
        'Enable','on',...
        'Tag','absolute',...
        'Position',[1 4.25 18 1.5]); 
hSD = uicontrol('parent',hth,'Style','Radio','Units','characters',...
        'String','Intensity std',...
        'Enable','on',...
        'Tag','SD',...
        'Position',[1 3 18 1.5]);
hAdaptive = uicontrol('parent',hth,'Style','Radio','Units','characters',...
        'String','Adaptive from ROIs',...
        'Enable','on',...
        'Tag','adaptive',...
        'Position',[1 1.75 18 1.5]);
hLegacySD = uicontrol('parent',hth,'Style','Radio','Units','characters',...
        'String','Intensity std (Legacy)',...
        'Enable','on',...
        'Tag','legacySD',...
        'Position',[1 0.5 18 1.5]);
    
if strcmp(thresh.units,'SD')
    set(hth,'SelectedObject',hSD);
elseif strcmp(thresh.units,'legacySD')
    set(hth,'SelectedObject',hLegacySD);
elseif strcmp(thresh.units,'adaptive')
    set(hth,'SelectedObject',hAdaptive);
else
    set(hth,'SelectedObject',hAbs);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%% contrast panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph = uipanel(fh,'Title','Contrast','Units','characters',...
             'Position',[30 46.2 83 8.5],'TitlePosition','centertop');

%info relative to the current slice        
u1 = uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','current z-plane',...
                'Position',[10,6.5,18,1.2]);        

u2 = uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min',...
                'Position',[0,5.4,8,1.2]);            
u3 = uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','median',...
                'Position',[13,5.4,8,1.2]);            
u4 = uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max',...
                'Position',[24,5.4,8,1.2]);            
u5 = uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','SD',...
                'Position',[39,5.4,8,1.2]);   
            
u6 = uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imin_Info_val',...
                'Style','text','String',strcurImin,...
                'Position',[0,4.6,8,1.2]);            
u7 =uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imed_Info_val',...
                'Style','text','String',strcurImed,...
                'Position',[13,4.6,8,1.2]);            
u8 = uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imax_Info_val',...
                'Style','text','String',strcurImax,...
                'Position',[24,4.6,8,1.2]);  
u9 = uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Isd_Info_val',...
                'Style','text','String',strcurIsd,...
                'Position',[39,4.6,8,1.2]);  

if nd == 2
    set(u1,'visible','off');
    set(u2,'visible','off');
    set(u3,'visible','off');
    set(u4,'visible','off');
    set(u5,'visible','off');
    set(u6,'visible','off');
    set(u7,'visible','off');
    set(u8,'visible','off');
    set(u9,'visible','off'); 
end

%info on whole stack                        
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','whole stack',...
                'Position',[10,2.3,18,1.2]);        

uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min',...
                'Position',[0,1.2,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','median',...
                'Position',[13,1.2,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max',...
                'Position',[24,1.2,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','SD',...
                'Position',[39,1.2,8,1.2]);
            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Imin_Info_val',...
                'Style','text','String',strImin,...
                'Position',[0,1,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Imed_Info_val',...
                'Style','text','String',strImed,...
                'Position',[13,1,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Imax_Info_val',...
                'Style','text','String',strImax,...
                'Position',[24,1,8,1.2]);                       
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Isd_Info_val',...
                'Style','text','String',strIsd,...
                'Position',[39,1,8,1.2]);             
            
%contrast settings           
u1 = uicontrol('Parent',ph,'Units','characters',...
                'Style','pushbutton',...
                'Tag','hslice_auto',...
                'String','auto adjust on current z-plane',...
                'Units','characters','Position',[52,3,30,1.9]); 
if nd == 2
    set(u1,'visible','off');
end

uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglob_auto',...
                'Style','pushbutton',...
                'String','auto adjust on whole stack',...
                'Units','characters','Position',[52,5,30,1.9]); 
            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min Int.',...
                'Position',[52,1.3,12,1.2]);
       
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hImin',...
                'Style','edit','String',strImin,...
                'Position',[52,0.4,12,1.2]);         

uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max Int.',...
                'Position',[66,1.3,12,1.2]);

uicontrol('Parent',ph,'Units','characters',...
                'Tag','hImax',...
                'Style','edit','String',strImax,...
                'Position',[66,0.4,12,1.2]);                
       
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% pointer info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uicontrol('Style','text','String','x',...
            'Units','characters',...
           'Position',[10,48.5,6,1.2]);

uicontrol('style','list','String',...
{'<HTML><FONT COLOR=00FF00>grab</FONT><HTML>',...
'<HTML><FONT COLOR=FF0000>snap</FONT><HTML>'},...
'Units','characters',...
'Tag','hstate',...
'Position',[120,53.3,12,1.2]);       

uicontrol('Style','text','String',strXpos,...
            'Units','characters',...
            'Tag','hxpos',...
           'Position',[14,48.5,12,1.2]);       
       
uicontrol('Style','text','String','y',...
           'Units','characters',...
           'Position',[10,47.3,6,1.2]);

uicontrol('Style','text','String',strYpos,...
            'Tag','hypos',...
            'Units','characters',...
           'Position',[14,47.3,12,1.2]);       
       
uicontrol('Style','text','String','Int',...
            'Units','characters',...
           'Position',[10,46.2,6,1.2]);

uicontrol('Style','text','String',strIval,...
            'Units','characters',...
            'Tag','hIval',...
           'Position',[14,46.2,12,1.2]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% z browser %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
u1 = uicontrol(fh,'Style','edit',...
                'Tag','hzpos',...
                'Units','characters',...
                'String',strZpos,...
                'Position',[54 3.5 6 1.2]);
if nd ==2
    zSliderStep = [1,1];
else
    zSliderStep = [1/(nz-1), 1/(nz-1)];
end
u2 = uicontrol(fh,'Style','slider',...
                'Tag','zslider',...
                'Units','characters',...
                'Max',nz,'Min',1,'Value',1,...
                'SliderStep',zSliderStep,...
                'Position',[10 2.3 102.4 1.2]);
            
if nd == 2
    set(u1,'visible','off');
    set(u2,'visible','off');
end

if nd == 2 || 3
    uicontrol(fh,'Style','text','Units','characters',...
        'String',['In grab mode, use arrows for fine (x,y) motion,',...
        'and the < or > keys for fine z motion'],...
        'HorizontalAlignment','center',...
        'Position',[52 0.4 52.4 1.2]);
else
    
    uicontrol(fh,'Style','text','Units','characters',...
        'String','In grab mode, use arrows for fine (x,y) motion',...
        'HorizontalAlignment','center',...
        'Position',[52 0.4 52.4 1.2]);
end

%% action button
uicontrol(fh,'Style','pushbutton',...
                    'Tag','hclose',...
                    'Units','characters',...
                    'String','done',...
                    'Position',[116 0.4 22 3]);
%% %%%%%%%%%%%%%%%%%%% overlay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uicontrol(fh,'Style','togglebutton',...
                    'Tag','hoverlay',...
                    'Units','characters',...
                    'Value',0,...
                    'String','Hide Overlay (shortcut: h)',...
                    'Position',[10 0.2 42 1.5]);                
                
                             
clear('nz','strXpos','strYpos','strZpos',...
    'Imed','Isd','Imin','Imax',...
    'strImin','strImax','strImed','strIsd',...
    'curImin','curImax','curImed','curIsd',...
    'strcurImin','strcurImax','strcurImed','strcurIsd',...
    'xpix','ypix','zpix',...
    'startzoomwidth','ph');

clear('stack1','stack1name');
end

