function varargout = explore_stack_v6(params,alData)


%% main figure
fh = figure('name',alData.curImgFile,...
              'Units','characters',...
              'MenuBar','none',...
              'Toolbar','none',...
              'NumberTitle','off',...
              'Position',[20 3 200 55],...
              'Visible','off'); 

%% parsing arguments 
if params.numDim == 3 || alData.isMovie
    [nx, ny, nz] = size(alData.img);
else
    [nx, ny] = size(alData.img);
    nz = 1;
end
setappdata(fh,'stackname',alData.curImgFile);
setappdata(fh,'nx',nx);
setappdata(fh,'ny',ny);
setappdata(fh,'nz',nz);
xfreeze = round(nx/2); 
yfreeze = round(ny/2); 
setappdata(fh,'xfreeze',xfreeze);
setappdata(fh,'yfreeze',yfreeze);
setappdata(fh,'fitres',[0,0,0,0,0,0,0]);
setappdata(fh,'defaultBackgroundColor','k');
setappdata(fh,'hIminVal',min(alData.img(:)));
setappdata(fh,'hImaxVal',max(alData.img(:)));
fh = set_figure_and_panel_for_explore_stack(fh,alData);
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
set(handles.hgauss_fit,'Callback',{@local_gaussian_fit_from_viewer,...
    alData,params,fh});
set(handles.hrecord,'Callback',{@record_fit_results,alData,fh});
set(handles.hxplot_range,'Callback',{@change_local_range,alData,fh});
set(handles.hyplot_range,'Callback',{@change_local_range,alData,fh});
set(handles.hzplot_range,'Callback',{@change_local_range,alData,fh});
set(fh,'WindowKeyPressFcn',{@moveCursorWithArrows,alData,fh});

set(fh,'SelectionType','alt');
set(fh,'CurrentAxes',handles.ha);
    
plot_current_z_stack_single_channel(fh,alData,1);
set(fh,'Visible','on');

uiwait;

varargout{1} = getappdata(fh,'fitres');
close(fh);
drawnow;
%% housekeeping
clear('nx','ny','nz','stack1','stack1name','xfreeze','yfreeze');
clear('fh','ha','handles');

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function moveCursorWithArrows(src,eventData,alData,fh)
    % list of keys that when pressed, map to actions
    keyList = {'uparrow','downarrow','leftarrow','rightarrow'};
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

        plot_local_data_single_channel(xfreeze,yfreeze,z,fh,alData);
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

function plot_current_z_stack_single_channel(fh,alData,z)

handles = guihandles(fh);    
%% parsing dimensions of input stack    
    if ndims(alData.img)==3
        [nx, ny, nz] = size(alData.img);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = alData.img(:,:,z);
        
    elseif ndims(alData.img) == 2
        [nx, ny] = size(alData.img);
        stack1 = alData.img;
    else
        disp('wrong stack size');
    end
    
    Imin = getCleanNumberFromEditEntry(handles,'hImin',fh,'hIminVal');
    Imax = getCleanNumberFromEditEntry(handles,'hImax',fh,'hImaxVal');

%% filling data 

    if(Imin >= Imax)
        stack1 = 0;
    else
        stack1 = min( max((stack1-Imin)/(Imax-Imin),0) , 1);
    end

%%    
    ha = handles.ha;
    set(fh,'CurrentAxes',handles.ha);
    
    imagesc(stack1); 
    colormap(gray);
    xlim([1,max(nx,ny)]);
    ylim([1,max(nx,ny)]);
    defaultBackgroundColor = getappdata(fh,'defaultBackgroundColor');
    set(ha,'Color',defaultBackgroundColor);
    set(ha,'Tag','ha');

clear('cur_slice','nx','ny','nz','Imin','Imax','newXLim','newYLim','stack1');
clear('stack1','handles','ha');
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

function plot_zoom_single_channel(x,y,z,width,alData,fh)
    
    handles = guihandles(fh); 
    zoomh = handles.zoomh;
    set(fh,'CurrentAxes',handles.zoomh);
    if ndims(alData.img) == 3
        [nx, ny, nz] = size(alData.img);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = alData.img(:,:,z); 
        
    elseif ndims(alData.img) ==2
        [nx, ny] = size(alData.img);
        stack1 = alData.img;
    end
    Imin = getCleanNumberFromEditEntry(handles,'hImin',fh,'hIminVal');
    Imax = getCleanNumberFromEditEntry(handles,'hImax',fh,'hImaxVal');
    
    %computing the value of the local slice 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(Imin >= Imax)
        stack1(1:nx,1:ny) = 0;
    else
        stack1 = min(1, max( (stack1-Imin) / (Imax-Imin),0));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %chosing the window corresponding to the cursor position
    if 2*width+1>nx
        if 2*width+1>ny
            imagesc(stack1);
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width y+width],[1 nx],stack1(1:nx,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[1,nx],stack1(1:nx,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[1,nx],stack1(1:nx,yy-width:yy+width));
        end
    elseif (x - width >= 1) && (x + width <= nx)
        if 2*width+1>ny
            imagesc([1,ny],[x-width,x+width],stack1(x-width:x+width,1:ny));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[x-width,x+width],stack1(x-width:x+width,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[x-width,x+width],stack1(x-width:x+width,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[x-width,x+width],stack1(x-width:x+width,yy-width:yy+width));
        end
    elseif (x - width < 1) && (x + width <= nx) 
        xx = width+1;
        if 2*width+1>ny
            imagesc([1,ny],[xx-width,xx+width],stack1(xx-width:xx+width,1:ny));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[xx-width,xx+width],stack1(xx-width:xx+width,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[xx-width,xx+width],stack1(xx-width:xx+width,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[xx-width,xx+width],stack1(xx-width:xx+width,yy-width:yy+width));
        end
    elseif (x + width > nx) && (x - width >= 1)
        xx = nx - width;
        if 2*width+1>ny
            imagesc([1,ny],[xx-width,xx+width],stack1(xx-width:xx+width,1:ny));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[xx-width,xx+width],stack1(xx-width:xx+width,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[xx-width,xx+width],stack1(xx-width:xx+width,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[xx-width,xx+width],stack1(xx-width:xx+width,yy-width:yy+width));
        end
    end 
    set(zoomh,'Tag','zoomh');
    colormap(gray);
    
    clear('stack1','nx','ny','nz','Imin','Imax','stack1','handles');
    clear('xx','yy');
    clear('x','y','z','width','zoomh','stack1','fh','hImin','hImax');
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

function change_local_range(src,eventdata,alData,fh)

    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_local_data_single_channel(x,y,z,fh,alData);
    clear('x','y','z','handles');
end

function record_fit_results(src,eventdata,alData,fh)
    handles = guihandles(fh);
    nd = ndims(alData.img);
    fit_res = getappdata(fh,'fitres');
    cur_res(1,1) = str2double(get(handles.hx0,'String'));
    cur_res(1,2) = str2double(get(handles.hy0,'String'));
    if (nd == 3) && ~alData.isMovie
        cur_res(1,3) = str2double(get(handles.hz0,'String'));
        cur_res(1,4) = str2double(get(handles.hI0,'String'));
        cur_res(1,5) = str2double(get(handles.hbg0,'String'));
        cur_res(1,6) = str2double(get(handles.hsxy0,'String'));
        cur_res(1,7) = str2double(get(handles.hsz0,'String'));
    else
        cur_res(1,3) = str2double(get(handles.hI0,'String'));
        cur_res(1,4) = str2double(get(handles.hbg0,'String'));
        cur_res(1,5) = str2double(get(handles.hsxy0,'String'));
    end
    
    n = str2double(get(handles.hFitsRecorded,'String'));
    if n == 0
        setappdata(fh,'fitres',cur_res);
    else
        setappdata(fh,'fitres',[fit_res;cur_res]);
    end
    set(handles.hFitsRecorded,'String',num2str(n+1));
    clear('handles','fit_res','cur_res','n');
end

function clean_close(src,eventdata,fh)
    uiresume;
end

function local_gaussian_fit_from_viewer(src,eventdata,alData,params,fh)

    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    handles = guihandles(fh);
    
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    
    nd = ndims(alData.img);
    if nd == 3
        Gaussout = gaussian_fit_local3(...
            alData,[xfreeze,yfreeze,z],params,1);
    else
        Gaussout = gaussian_fit_local3(...
            alData,[xfreeze,yfreeze],params,1);
    end
    if nd == 3 && ~alData.isMovie
        set(handles.hI0,'String',num2str(Gaussout(1)));
        set(handles.hbg0,'String',num2str(Gaussout(2)));
        set(handles.hsxy0,'String',num2str(Gaussout(3)));
        set(handles.hsz0,'String',num2str(Gaussout(4)));
        set(handles.hx0,'String',num2str(Gaussout(5)));
        set(handles.hy0,'String',num2str(Gaussout(6)));
        set(handles.hz0,'String',num2str(Gaussout(7)));
    else
        set(handles.hI0,'String',num2str(Gaussout(1)));
        set(handles.hbg0,'String',num2str(Gaussout(2)));
        set(handles.hsxy0,'String',num2str(Gaussout(3)));
        set(handles.hx0,'String',num2str(Gaussout(4)));
        set(handles.hy0,'String',num2str(Gaussout(5)));
    end
    
    plot_local_data_and_fit_single_channel_bgcorr(...
        xfreeze,yfreeze,z,fh,alData,params,Gaussout);
    clear('xfreeze','yfreeze','handles','z','p','Gaussout','resnorm');
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
    
    val = get(handles.hstate,'Value');
    val = 3-val; %1->2 or 2->1
    set(handles.hstate,'Value',val);
    
    if val == 2   %if switching to snap state
        set(fh,'CurrentAxes',handles.ha);
        point = get(gca,'CurrentPoint');
        yfreeze = round(point(1,1));    %note the inverse convention for mouse and array (xy -> yx)
        xfreeze = round(point(1,2));
        xfreeze = max(1,min(xfreeze,nx));
        yfreeze = max(1,min(yfreeze,ny));
        
    end
    setappdata(fh,'xfreeze',xfreeze);
    setappdata(fh,'yfreeze',yfreeze);
    
    clear('point','val','point','xfreeze','yfreeze');
    clear('nx','ny','handles');
end

function slider_change(src,eventdata,alData,fh)
    
    handles = guihandles(fh);
    
    %updating the z-slider and the z info window
    z = get(handles.zslider,'Value');
    z = round(z);
    set(handles.zslider,'Value',z);
    set(handles.hzpos,'String',num2str(z));
    
    %updating the Imin/Imax info
    [curImin,curImax,curImed,curIsd] = getIntensityStats(alData,z);
    set(handles.hslice_Imin_Info_val,'String',num2str(curImin));
    set(handles.hslice_Imed_Info_val,'String',num2str(curImed));
    set(handles.hslice_Imax_Info_val,'String',num2str(curImax)); 
    set(handles.hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    Int = alData.img(x,y,z);
    set(handles.hIval,'String',num2str(Int));

    plot_local_data_single_channel(x,y,z,fh,alData);
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,alData,z);  
    
    clear('x','y','z','Int','curImin','curImax','curImed','curIsd');
    clear('nx','ny','xfreeze','yfreeze','handles');
end

function zpos_change(src,eventdata,alData,fh)

    nz = getappdata(fh,'nz');
    
    handles = guihandles(fh);
    
    %updating the z-slider and the z info window
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    z = max(min(z,nz) , 1);
    set(handles.hzpos,'String',num2str(z));
    set(handles.zslider,'Value',z);
    
    %updating the Imin/Imax info
    [curImin,curImax,curImed,curIsd] = getIntensityStats(alData,z);
    set(handles.hslice_Imin_Info_val,'String',num2str(curImin));
    set(handles.hslice_Imed_Info_val,'String',num2str(curImed));
    set(handles.hslice_Imax_Info_val,'String',num2str(curImax));
    set(handles.hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,alData,z);
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    Int = alData.img(x,y,z);
    set(handles.hIval,'String',num2str(Int));

    plot_local_data_single_channel(x,y,z,fh,alData);
    
    clear('x','y','z','Int','curImin','curImax','curImed','curIsd'); 
    clear('nx','ny','xfreeze','yfreeze','handles');
end

function [Imin,Imax,Imed,Isd] = getIntensityStats(alData,z)
    nd = ndims(alData.img);
    if nd == 3
        Imin = min(alData.img(:,:,z),[],'all'); 
        Imax = max(alData.img(:,:,z),[],'all'); 
        Imed = median(alData.img(:,:,z),'all'); 
        Isd = std(alData.img(:,:,z),0,'all');
    else
        Imin = min(alData.img(:,:),[],'all'); 
        Imax = max(alData.img(:,:),[],'all'); 
        Imed = median(alData.img(:,:),'all'); 
        Isd = std(alData.img(:,:),0,'all');
    end
end

function change_Iminmax(src,eventdata,alData,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,alData,z);
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_zoom_single_channel(x,y,z,zoomwidth,alData,fh);
    
    clear('x','y','z','zoomwidth');
    clear('handles');
end

function slice_auto_adjust(src,eventdata,alData,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %updating Imin / Imax info
    [Imin,Imax,~,~] = getIntensityStats(alData,z);
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,alData,z);
    
    %replotting the zoom window
    zoomwidth = str2douvble(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_single_channel(x,y,z,zoomwidth,alData,fh);
    clear('x','y','z','zoomwidth','Imin','Imax');
    clear('handles');
end

function global_auto_adjust(src,eventdata,alData,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %updating the Imin/Imax info
    [Imin,Imax,~,~] = getIntensityStats(alData,z);
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,alData,z);
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_single_channel(x,y,z,zoomwidth,alData,fh);
    
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
    nz = getappdata(fh,'nz');
    
    handles = guihandles(fh);
    val = get(handles.hstate,'Value');  %this is the snap/grab mode
   
    %I replot things only if the grab mode is active
    if val == 1 

        %I get the current (x,y,z) 
        [x,y] = get_current_pointer_position(fh); 
      
        z = get(handles.zslider,'Value');
        z = round(z);
        
        x = max(1,min(x,nx));
        y = max(1,min(y,ny));
        z = max(1,min(z,nz));
        
        %I update the current display of the pointer info
        set(handles.hxpos,'String',num2str(x));
        set(handles.hypos,'String',num2str(y));
        Int = alData.img(x,y,z);
        set(handles.hIval,'String',num2str(Int));
        
        %plotting the profiles and zoom window
        plot_local_data_single_channel(x,y,z,fh,alData);
        
    end
    clear('x','y','z','val','Int','nx','ny','handles');
    
    % this command is key to ensure that there is no pile up of callbacks
    % which leads to lags
    drawnow;
end

function plot_local_data_single_channel(x,y,z,fh,alData)
    
    handles = guihandles(fh);
    
    xrange = str2double(get(handles.hxplot_range,'String'));
    yrange = str2double(get(handles.hyplot_range,'String'));
    nd = ndims(alData.img);
    if nd == 3
        zrange = str2double(get(handles.hzplot_range,'String'));
    end
    zoomwidth = str2double(get(handles.hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile(x,y,z,xrange,alData,fh);
    plot_yprofile(x,y,z,yrange,alData,fh);
    if nd == 3
        plot_zprofile(x,y,z,zrange,alData,fh);
    end
    plot_zoom_single_channel(x,y,z,zoomwidth,alData,fh);
    
    clear('handles','xrange','yrange','zrange','zoomwidth');
end

function plot_local_data_and_fit_single_channel_bgcorr(...
    x,y,z,fh,alData,params,Gaussout)
    
    handles = guihandles(fh);
    xrange = str2double(get(handles.hxplot_range,'String'));
    yrange = str2double(get(handles.hyplot_range,'String'));
    nd = ndims(alData.img);
    if nd == 3
        zrange = str2double(get(handles.hzplot_range,'String'));
    end
    zoomwidth = str2double(get(handles.hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile_and_fit_bgcorr(x,y,z,xrange,alData,fh,params,Gaussout);
    plot_yprofile_and_fit_bgcorr(x,y,z,yrange,alData,fh,params,Gaussout);
    if nd == 3
        plot_zprofile_and_fit_bgcorr(x,y,z,zrange,alData,fh,params,Gaussout);
    end
    plot_zoom_single_channel(x,y,z,zoomwidth,alData,fh);
    
    clear('xrange','yrange','zrange','zoomwidth','handles');
    
end

function plot_xprofile_and_fit_bgcorr(xarr,yarr,z,xrange,alData,fh,params,Gaussout)
    handles = guihandles(fh);
    xh = handles.xh;
    nx = size(alData.img,1);
    nd = ndims(alData.img);
    if nd == 3
        xprofile = alData.img(1:nx,yarr,z);
    else
        xprofile = alData.img(1:nx,yarr);
    end
    switch params.psfType
        case 'gaussian'
            if nd == 3 && ~alData.isMovie
                xfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_gaussian3D((1:nx)',yarr,z,...
                    Gaussout(5),Gaussout(6),Gaussout(7),Gaussout(3),...
                    Gaussout(4),1,1,1);
            else
                xfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_gaussian2D((1:nx)',yarr,...
                    Gaussout(4),Gaussout(5),Gaussout(3),1,1);
            end
        case 'integratedGaussian'
            if nd == 3 && ~alData.isMovie
                I0 = intensity_integrated_gaussian3D(...
                    1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
                xfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_integrated_gaussian3D((1:nx)',...
                    yarr,z,Gaussout(5),Gaussout(6),Gaussout(7),...
                    Gaussout(3),Gaussout(4),1,1,1)/I0;
            else
                I0 = intensity_integrated_gaussian2D(...
                1,1,0.5,0.5,Gaussout(3),1,1);
                xfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_integrated_gaussian2D((1:nx)',...
                    yarr,Gaussout(4),Gaussout(5),Gaussout(3),1,1)/I0;
            end
        case 'integratedGaussianStdZ'
            if nd == 3 && ~alData.isMovie
                I0 = intensity_integrated_gaussian3D_stdZ(...
                    1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
                xfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian3D_stdZ((1:nx)',yarr,z,...
                 Gaussout(5),Gaussout(6),Gaussout(7),Gaussout(3),Gaussout(4),1,1,1)/I0;
            elseif nd == 3 && alData.isMovie
                I0 = intensity_integrated_gaussian2D(...
                    1,1,0.5,0.5,Gaussout(3),1,1);
                xfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_integrated_gaussian2D((1:nx)',...
                    yarr,Gaussout(4),Gaussout(5),Gaussout(3),1,1)/I0;
            else
                disp('Cannot use integratedGaussianStdZ psfType on 2D data');
            end
    end
    
    % add background
    if nd == 3 && ~alData.isMovie
        xfit = xfit + ...
            Gaussout(10)*(1:nx)' + ...
            Gaussout(11)*yarr + Gaussout(12)*z + Gaussout(13);
    else
        xfit = xfit + ...
            Gaussout(8)*(1:nx)' + Gaussout(9)*yarr + Gaussout(10);
    end
    
    x1 = max(1,xarr-xrange);
    x2 = min(nx,xarr+xrange);
    Imin = min(xprofile(x1:x2));
    Imax = max(xprofile(x1:x2));
    
    set(fh,'CurrentAxes',handles.xh);
    plot(x1:x2,xprofile(x1:x2),[xarr xarr],[Imin Imax]); 
    set(gca,'nextPlot','add');
    plot(x1:x2,xfit(x1:x2),'r',[xarr xarr],[Imin Imax]); 
    set(gca,'nextPlot','replace');
    set(xh,'Tag','xh');
    clear('nx','xprofile','Imin','Imax','x','pos','xfit');
    clear('handles');
end

function plot_yprofile_and_fit_bgcorr(xarr,yarr,z,yrange,alData,fh,params,Gaussout)
    
    handles = guihandles(fh);
    yh = handles.yh;
    ny = size(alData.img,2);
    nd = ndims(alData.img);
    if nd == 3 
        yprofile = alData.img(xarr,1:ny,z);
    else
        yprofile = alData.img(xarr,1:ny);
    end

    switch params.psfType
        case 'gaussian'
            if nd == 3 && ~alData.isMovie
                yfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_gaussian3D(xarr,(1:ny)',z,...
                    Gaussout(5),Gaussout(6),Gaussout(7),Gaussout(3),...
                    Gaussout(4),1,1,1);
            else
                yfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_gaussian2D(xarr,(1:ny)',...
                    Gaussout(4),Gaussout(5),Gaussout(3),1,1);
            end
        case 'integratedGaussian'
            if nd ==3 && ~alData.isMovie
                I0 = intensity_integrated_gaussian3D(...
                    1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
                yfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_integrated_gaussian3D(...
                    xarr,(1:ny)',z,Gaussout(5),Gaussout(6),Gaussout(7),...
                    Gaussout(3),Gaussout(4),1,1,1)/I0;
            else
                I0 = intensity_integrated_gaussian2D(...
                    1,1,0.5,0.5,Gaussout(3),1,1);
                yfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian2D(xarr,(1:ny)',...
                    Gaussout(4),Gaussout(5),Gaussout(3),1,1)/I0;
            end
        case 'integratedGaussianStdZ'
            if nd ==3 && ~alData.isMovie
                I0 = intensity_integrated_gaussian3D_stdZ(...
                    1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
                yfit = Gaussout(2)+ ...
                    Gaussout(1)*intensity_integrated_gaussian3D_stdZ(...
                    xarr,(1:ny)',z,Gaussout(5),Gaussout(6),Gaussout(7),...
                    Gaussout(3),Gaussout(4),1,1,1)/I0;
             else
                I0 = intensity_integrated_gaussian2D(...
                    1,1,0.5,0.5,Gaussout(3),1,1);
                yfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian2D(xarr,(1:ny)',...
                    Gaussout(4),Gaussout(5),Gaussout(3),1,1)/I0;
            end
    end
    if nd == 3 && ~alData.isMovie
        yfit = yfit + ...
            Gaussout(10)*xarr + ...
            Gaussout(11)*(1:ny)' + Gaussout(12)*z + Gaussout(13);
    else
        yfit = yfit + ...
            Gaussout(8)*xarr + Gaussout(9)*(1:ny)' + Gaussout(10);
    end
    y1 = max(1,yarr-yrange);
    y2 = min(ny,yarr+yrange);
    Imin = min(yprofile(y1:y2));
    Imax = max(yprofile(y1:y2));
    
    set(fh,'CurrentAxes',handles.yh);
    plot(y1:y2,yprofile(y1:y2),[yarr yarr],[Imin Imax]); 
    %hold(handles.yh,'on');
    set(gca,'nextPlot','add');
    plot(y1:y2,yfit(y1:y2),'r',[yarr yarr],[Imin Imax]); 
    %hold(handles.yh,'off'); 
    set(gca,'nextPlot','replace');
    
    set(yh,'Tag','yh');
    clear('ny','yprofile','Imin','Imax','y','pos','yfit');
    clear('xarr','yarr','handles');
end

function plot_zprofile_and_fit_bgcorr(xarr,yarr,z,zrange,alData,fh,params,Gaussout)
    
    nd = ndims(alData.img);
    if nd ~= 3
        return
    end
    handles = guihandles(fh);
    zh = handles.zh;
    nz = size(alData.img,3);
    zprofile = reshape(alData.img(xarr,yarr,1:nz),1,nz);
    if alData.isMovie
        z1 = max(1,z-zrange);
        z2 = min(nz,z+zrange);
        Imin = min(zprofile(z1:z2));
        Imax = max(zprofile(z1:z2));
        set(fh,'CurrentAxes',handles.zh);
        plot(z1:z2,zprofile(z1:z2),[z z],[Imin Imax]); 
        set(gca,'nextPlot','replace');
        set(zh,'Tag','zh');
        return
    end
    
    switch params.psfType
        case 'gaussian'
            zfit = Gaussout(2)+ Gaussout(1)*intensity_gaussian3D(xarr,yarr,(1:nz)',...
                Gaussout(5),Gaussout(6),Gaussout(7),Gaussout(3),Gaussout(4),1,1,1);
        case 'integratedGaussian'
            I0 = intensity_integrated_gaussian3D(...
                1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
            zfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian3D(xarr,yarr,(1:nz)',...
                Gaussout(5),Gaussout(6),Gaussout(7),Gaussout(3),Gaussout(4),1,1,1)/I0;
        case 'integratedGaussianStdZ'
            I0 = intensity_integrated_gaussian3D_stdZ(...
                1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
            zfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian3D_stdZ(xarr,yarr,(1:nz)',...
                Gaussout(5),Gaussout(6),Gaussout(7),Gaussout(3),Gaussout(4),1,1,1)/I0;
    end
    
    zfit = zfit + Gaussout(10)*xarr + ...
            Gaussout(11)*yarr + Gaussout(12)*(1:nz)' + Gaussout(13);
    z1 = max(1,z-zrange);
    z2 = min(nz,z+zrange);
    Imin = min(zprofile(z1:z2));
    Imax = max(zprofile(z1:z2));
    
    set(fh,'CurrentAxes',handles.zh);
    plot(z1:z2,zprofile(z1:z2),[z z],[Imin Imax]); 
    %hold(handles.zh,'on');
    set(gca,'nextPlot','add');
    plot(z1:z2,zfit(z1:z2),'r',[z z],[Imin Imax]); 
    %hold(handles.zh,'off');
    set(gca,'nextPlot','replace');
    
    set(zh,'Tag','zh');
    clear('nz','zprofile','Imin','Imax','z','zz','pos','zfit');
    clear('xarr','yarr','handles');
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
    nx = size(alData.img,1);
    nd = ndims(alData.img);
    if nd == 3
        xprofile = alData.img(1:nx,yarr,z);
    else
        xprofile = alData.img(1:nx,yarr);
    end
    xh = handles.xh;
    x1 = max(1,xarr-xrange);
    x2 = min(nx,xarr+xrange);
    Imin = min(xprofile(x1:x2));
    Imax = max(xprofile(x1:x2));
    
    set(fh,'CurrentAxes',handles.xh);
    plot(x1:x2,xprofile(x1:x2),[xarr xarr],[Imin Imax]); 
    set(xh,'Tag','xh');
    clear('nx','xprofile','Imin','Imax','x','handles','xh');
end

function plot_yprofile(xarr,yarr,z,yrange,alData,fh)

    handles = guihandles(fh);
    ny = size(alData.img,2);
    nd = ndims(alData.img);
    if nd == 3
        yprofile = alData.img(xarr,1:ny,z);
    else
        yprofile = alData.img(xarr,1:ny);
    end
    yh = handles.yh;
    y1 = max(1,yarr-yrange);
    y2 = min(ny,yarr+yrange);
    Imin = min(yprofile(y1:y2));
    Imax = max(yprofile(y1:y2));
    
    set(fh,'CurrentAxes',handles.yh);
    plot(y1:y2,yprofile(y1:y2),[yarr yarr],[Imin Imax]); 
    set(yh,'Tag','yh');
    clear('ny','yprofile','Imin','Imax','handles','yh');
end

function plot_zprofile(xarr,yarr,z,zrange,alData,fh)
    nd = ndims(alData.img);
    if nd ~= 3
        return
    end
    handles = guihandles(fh);
    nz = size(alData.img,3); 
    zprofile = reshape(alData.img(xarr,yarr,1:nz),1,nz);
    
    zh = handles.zh;
    z1 = max(1,z-zrange);
    z2 = min(nz,z+zrange);
    Imin = min(zprofile(z1:z2));
    Imax = max(zprofile(z1:z2));
    
    set(fh,'CurrentAxes',handles.zh);
    plot(z1:z2,zprofile(z1:z2),[z z],[Imin Imax]); 
    set(zh,'Tag','zh');
    clear('nz','zprofile','Imin','Imax','handles','zh');
end

function fh = set_figure_and_panel_for_explore_stack(fh,alData)


%%%%%%%%%%%%%%%%%%%%% variables initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz = getappdata(fh,'nz');
nd = ndims(alData.img);

strXpos = '1';     % default x value of the current pixel
strYpos = '1';     % default y value of the current pixel
if nd == 3
    strIval = num2str(alData.img(1,1,1));   % value of the intensity @ the current pixel
else
    strIval = num2str(alData.img(1,1)); 
end
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

%values of min max and median of current slice (z=1 when you open the viewer)
[curImin,curImax,curImed,curIsd] = getIntensityStats(alData,1);

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

set(fh,'Name',getappdata(fh,'stackname'));              
axes('Parent',fh,'Units','characters','Color', 'k',...
    'Position',[10,6.15,102.4,39.4],'Tag','ha');     %main figure
        
axes('Parent',fh,'Units','characters','Position',[120 45.4 66 7.7],'Tag','xh'); 
axes('Parent',fh,'Units','characters','Position',[120 33.8 66 7.7],'Tag','yh'); 
axes_zh = axes('Parent',fh,'Units','characters','Position',[120 22.3 66 7.7],'Tag','zh'); 
if nd == 2
    set(axes_zh,'visible','off');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% settings for the axial views %%%%%%%%%%%%%%%%
uicontrol('Units','characters',...
                'Style','text','String','profile along x',...
                'Position',[144,53.3,20,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','profile along y',...
                'Position',[144,41.8,20,1.2]);
uiZprofile = uicontrol('Units','characters',...
                'Style','text','String','profile along z',...
                'Position',[144,30.3,20,1.2]);
if nd == 2
    set(uiZprofile,'visible','off');
end            
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,49.2,10,3.6]); 
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,37.7,10,3.6]); 
uiZrange = uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,26.2,10,3.6]); 
if nd == 2
    set(uiZrange,'visible','off');
end
uicontrol('Units','characters',...
                'Tag','hxplot_range',...
                'Style','edit','String',xpix,...
                'Position',[190,47.7,6,1.2]);                  
uicontrol('Units','characters',...
                'Tag','hyplot_range',...
                'Style','edit','String',ypix,...
                'Position',[190,36.2,6,1.2]);            
uiZplotRange = uicontrol('Units','characters',...
                'Tag','hzplot_range',...
                'Style','edit','String',zpix,...
                'Position',[190,24.6,6,1.2]);
if nd == 2
    set(uiZplotRange,'visible','off');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% zoom plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
axes('Parent',fh,'Units','characters','Position',[144 1.9 46 17.7],'Tag','zoomh');             

uicontrol('Units','characters',...
                'Style','text','String','window size',...
                'Position',[192,12,8,2.3]);
            
uicontrol('Units','characters',...
                'Tag','hzoom_width',...
                'Style','edit','String',startzoomwidth,...
                'Position',[192,10.4,6,1.2]); 

%% %%%%%%%%%%%%%%%%%%%%%%%% fit gaussian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uicontrol('Units','characters',...
                'Tag','hgauss_fit',...
                'Style','pushbutton','String','local gaussian fit',...
                'Position',[116,16.7,22,2]);

uicontrol('Units','characters',...
                'Style','text','String','x max',...
                'HorizontalAlignment','right',...
                'Position',[116,15.4,8,1.2]);
uicontrol('Units','characters',...
                'Tag','hx0',...
                'Style','text','String','undefined',...
                'Position',[126,15.4,12,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','y max',...
                'HorizontalAlignment','right',...
                'Position',[116,14.2,8,1.2]);
uicontrol('Units','characters',...
                'Tag','hy0',...
                'Style','text','String','undefined',...
                'Position',[132,14.2,6,1.2]);            
uiZmax = uicontrol('Units','characters',...
                'Style','text','String','z max',...
                'HorizontalAlignment','right',...
                'Position',[116,13.1,8,1.2]);            
uihz0 = uicontrol('Units','characters',...
                'Tag','hz0',...
                'Style','text','String','undefined',...
                'Position',[132,13.1,6,1.2]);
if nd == 2
    set(uiZmax,'visible','off');
    set(uihz0,'visible','off');
end
uicontrol('Units','characters',...
                'Style','text','String','integrated Intensity',...
                'HorizontalAlignment','right',...
                'Position',[113,11.5,17,1.2]);
uicontrol('Units','characters',...
                'Tag','hI0',...
                'Style','text','String','undefined',...
                'Position',[132,11.5,6,1.2]); 
uicontrol('Units','characters',...
                'Style','text','String','background',...
                'HorizontalAlignment','right',...
                'Position',[113,10.4,17,1.2]);
uicontrol('Units','characters',...
                'Tag','hbg0',...
                'Style','text','String','undefined',...
                'Position',[132,10.4,6,1.2]); 
uicontrol('Units','characters',...
                'Style','text','String','s_xy',...
                'HorizontalAlignment','right',...
                'Position',[113,8.8,17,1.2]);
uicontrol('Units','characters',...
                 'Tag','hsxy0',...
                'Style','text','String','undefined',...
                'Position',[132,8.8,6,1.2]);  
uiSz = uicontrol('Units','characters',...
                'Style','text','String','s_z',...
                'HorizontalAlignment','right',...
                'Position',[113,7.7,17,1.2]);
uihsz0 = uicontrol('Units','characters',...
                'Tag','hsz0',...
                'Style','text','String','undefined',...
                'Position',[132,7.7,6,1.2]);             
if nd == 2
    set(uiSz,'visible','off');
    set(uihsz0,'visible','off');
end            
uicontrol('Units','characters',...
                'Tag','hrecord',...
                'Style','PushButton','String','Record Fit Results',...
                'Position',[116,5.5,22,2]);                       
uicontrol('Units','characters',...
                'Style','text','String','Fits Recorded',...
                'HorizontalAlignment','left',...
                'Position',[116,4.1,14,1.2]);
uicontrol('Units','characters',...
                'Tag','hFitsRecorded',...
                'Style','text','String','0',...
                'Position',[132,4.1,6,1.2]);            
                 
%% %%%%%%%%%%%%%%%%%%%%%%%%% contrast panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph = uipanel(fh,'Title','Contrast','Units','characters',...
             'Position',[30 46.2 83 8.5],'TitlePosition','centertop');

%info relative to the current slice        
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','current z-plane',...
                'Position',[10,6.5,18,1.2]);        

uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min',...
                'Position',[0,5.4,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','median',...
                'Position',[13,5.4,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max',...
                'Position',[24,5.4,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','SD',...
                'Position',[39,5.4,8,1.2]);   
            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imin_Info_val',...
                'Style','text','String',strcurImin,...
                'Position',[0,4.6,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imed_Info_val',...
                'Style','text','String',strcurImed,...
                'Position',[13,4.6,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imax_Info_val',...
                'Style','text','String',strcurImax,...
                'Position',[24,4.6,8,1.2]);  
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Isd_Info_val',...
                'Style','text','String',strcurIsd,...
                'Position',[39,4.6,8,1.2]);  
            
%info on whole stack 
if nd == 3
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
end   

%contrast settings           
uicontrol('Parent',ph,'Units','characters',...
                'Style','pushbutton',...
                'Tag','hslice_auto',...
                'String','auto adjust on current z-plane',...
                'Units','characters','Position',[52,3,30,1.9]); 

uiGlob = uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglob_auto',...
                'Style','pushbutton',...
                'String','auto adjust on whole stack',...
                'Units','characters','Position',[52,5,30,1.9]); 
if nd ==2
    set(uiGlob,'visible','off');
end
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
uiZpos = uicontrol(fh,'Style','edit',...
                'Tag','hzpos',...
                'Units','characters',...
                'String',strZpos,...
                'Position',[54 3.5 6 1.2]);
if nd ==2
    set(uiZpos,'visible','off');
end

if nd ==2
    zSliderStep = [1,1];
else
    zSliderStep = [1/(nz-1), 1/(nz-1)];
end
uiZslider = uicontrol(fh,'Style','slider',...
                'Tag','zslider',...
                'Units','characters',...
                'Max',nz,'Min',1,'Value',1,...
                'SliderStep',zSliderStep,...
                'Position',[10 2.3 102.4 1.2]);
if nd ==2
    set(uiZslider,'visible','off');
end

if nd ==2
    uicontrol(fh,'Style','text','Units','characters',...
        'String','In grab mode, use arrows for fine (x,y) motion',...
        'HorizontalAlignment','center',...
        'Position',[10 0.8 102.4 1.2]);
else
    uicontrol(fh,'Style','text','Units','characters',...
        'String',['In grab mode, use arrows for fine (x,y) motion,',...
        'and the < or > keys for fine z motion'],...
        'HorizontalAlignment','center',...
        'Position',[10 0.8 102.4 1.2]);
end

uicontrol(fh,'Style','pushbutton',...
                    'Tag','hclose',...
                    'Units','characters',...
                    'String','done',...
                    'Position',[116 0.4 22 3]);

clear('nz','strXpos','strYpos','strZpos',...
    'Imed','Isd','Imin','Imax',...
    'strImin','strImax','strImed','strIsd',...
    'curImin','curImax','curImed','curIsd',...
    'strcurImin','strcurImax','strcurImed','strcurIsd',...
    'xpix','ypix','zpix',...
    'startzoomwidth','ph');

end