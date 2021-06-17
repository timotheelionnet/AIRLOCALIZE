function [stack, img_read, filename] = tiffread5(filename, img_first, img_last)
% tiffread, version 2.4
%
% [stack, nbImages] = tiffread;
% [stack, nbImages] = tiffread(filename);
% [stack, nbImages] = tiffread(filename, imageIndex);
% [stack, nbImages] = tiffread(filename, firstImageIndex, lastImageIndex);
%
% Reads 8,16,32 bits uncompressed grayscale and (some) color tiff files,
% as well as stacks or multiple tiff images, for example those produced
% by metamorph or NIH-image. However, the entire TIFF standard is not
% supported (but you may extend it).
%
% The function can be called with a file name in the current directory,
% or without argument, in which case it pop up a file openning dialog
% to allow manual selection of the file.
% If the stacks contains multiples images, loading can be restricted by
% specifying the first and last images to read, or just one image to read.
%
% at return, nbimages contains the number of images read, and S is a vector
% containing the different images with some additional informations. The
% image pixels values are stored in the field .data, for gray level images,
% or in the fields .red, .green and .blue
% the pixels values are in the native (integer) format,
% and must be converted to be used in most matlab functions.
%
% Example:
% im = tiffread('spindle.stk');
% imshow( double(im(5).data), [] );
%
% Francois Nedelec, EMBL, Copyright 1999-2006.
% rewriten July 7th, 2004 at Woods Hole during the physiology course.
% last modified April 12, 2006.
% Contributions:
%   Kendra Burbank suggested the waitbar
%   Hidenao Iwai for the code to read floating point images,
%   Stephen Lang made tiffread more compliant with PlanarConfiguration
%
% Please, help us improve this software: send us feedback/bugs/suggestions
% This software is provided at no cost by a public research institution.
% However, postcards are always welcome!
%
% Francois Nedelec
% nedelec (at) embl.de
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
% http://www.embl.org
% http://www.cytosim.org

%Optimization: join adjacent TIF strips: this results in faster reads
consolidateStrips = 1;

%if there is no argument, we ask the user to choose a file:
if (nargin == 0)
    [filename, pathname] = uigetfile('*.tif;*.stk;*.lsm', 'select image file');
    filename = [ pathname, filename ];
end

if (nargin<=1);  img_first = 1; img_last = 10000; end
if (nargin==2);  img_last = img_first;            end


% not all valid tiff tags have been included, as they are really a lot...
% if needed, tags can easily be added to this code
% See the official list of tags:
% http://partners.adobe.com/asn/developer/pdfs/tn/TIFF6.pdf
%
% the structure IMG is returned to the user, while TIF is not.
% so tags usefull to the user should be stored as fields in IMG, while
% those used only internally can be stored in TIF.

global TIF;
TIF = [];

%% initializing

%counters for the number of images read and skipped
img_skip  = 0;
img_read  = 0;

% set defaults values :
TIF.SampleFormat     = 1;
TIF.SamplesPerPixel  = 1;
TIF.BOS              = 'l';          %byte order string

if  isempty(findstr(filename,'.'))
    filename = [filename,'.tif'];
end

TIF.file = fopen(filename,'r','l');
if TIF.file == -1
    filename = strrep(filename, '.tif', '.stk');
    TIF.file = fopen(filename,'r','l');
    if TIF.file == -1
        error(['file <',filename,'> not found.']);
    end
end


%% read header

%% read byte order: II = little endian, MM = big endian
byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
    TIF.BOS = 'l';                                %normal PC format
elseif ( strcmp(byte_order','MM') )
    TIF.BOS = 'b';
else
    error('This is not a TIFF file (no MM or II).');
end

%% read in a number which identifies file as TIFF format
tiff_id = fread(TIF.file,1,'uint16', TIF.BOS);
if (tiff_id ~= 42)
    error('This is not a TIFF file (missing 42).');
end

%% read the byte offset for the first image file directory (IFD)
ifd_pos = fread(TIF.file,1,'uint32', TIF.BOS);



while (ifd_pos ~= 0)
    
    clear IMG;
    IMG.filename = fullfile( pwd, filename );
    % move in the file to the first IFD
    fseek(TIF.file, ifd_pos, -1);
    %disp(strcat('reading img at pos :',num2str(ifd_pos)));

    
    %read in the number of IFD entries
    num_entries = fread(TIF.file,1,'uint16', TIF.BOS);
    %disp(strcat('num_entries =', num2str(num_entries)));

    %% read and process each IFD entry
    for i = 1:num_entries

        % save the current position in the file
        file_pos  = ftell(TIF.file);

        % read entry tag
        TIF.entry_tag = fread(TIF.file, 1, 'uint16', TIF.BOS);
        entry = readIFDentry;
        %disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));

        switch TIF.entry_tag
            case 256         % image width - number of column
                IMG.width          = entry.val;
            case 257         % image height - number of row
                IMG.height         = entry.val;
                TIF.ImageLength    = entry.val;
            case 258         % BitsPerSample per sample
                TIF.BitsPerSample  = entry.val;
                TIF.BytesPerSample = TIF.BitsPerSample / 8;
                IMG.bits           = TIF.BitsPerSample(1);
                %fprintf(1,'BitsPerSample %i %i %i\n', entry.val);
            case 273         % strip offset
                TIF.StripOffsets   = entry.val;
                TIF.StripNumber    = entry.cnt;
                %fprintf(1,'StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
            case 277         % sample_per pixel
                TIF.SamplesPerPixel  = entry.val;
                %fprintf(1,'Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
            case 278         % rows per strip
                TIF.RowsPerStrip   = entry.val;
                %fprintf(1,'RowsPerStrip = %i %i\n', TIF.RowsPerStrip);
            case 279         % strip byte counts - number of bytes in each strip after any compressio
                TIF.StripByteCounts= entry.val;
            case 284         %planar configuration describe the order of RGB
                TIF.PlanarConfiguration = entry.val;
            case 339
                TIF.SampleFormat   = entry.val;
            case 33629       %this tag identify the image as a Metamorph stack!
                TIF.MM_stack       = entry.val;
                TIF.MM_stackCnt    = entry.cnt;
                if ( img_last > img_first )
                    waitbar_handle = waitbar(0,'Please wait...','Name',['Reading ' filename]);
                end
            otherwise
                %fprintf(1,'ignored TIFF entry with tag %i (cnt %i)\n', TIF.entry_tag, entry.cnt);
        end
        % move to next IFD entry in the file
        fseek(TIF.file, file_pos+12,-1);
    end

    %Planar configuration is not fully supported
    %Per tiff spec 6.0 PlanarConfiguration irrelevent if SamplesPerPixel==1
    %Contributed by Stephen Lang
    if ((TIF.SamplesPerPixel ~= 1) && (TIF.PlanarConfiguration == 1))
        error('PlanarConfiguration = %i not supported', TIF.PlanarConfiguration);
    end

    %% total number of bytes per image:
    PlaneBytesCnt = IMG.width * IMG.height * TIF.BytesPerSample;

    %% Try to consolidate the strips into a single one to speed-up reading:  
    if consolidateStrips
        BytesCnt = TIF.StripByteCounts(1);

        if BytesCnt < PlaneBytesCnt
            ConsolidateCnt = 1;
            %Count how many Strip are needed to produce a plane
            while TIF.StripOffsets(1) + BytesCnt == TIF.StripOffsets(ConsolidateCnt+1)
                ConsolidateCnt = ConsolidateCnt + 1;
                BytesCnt = BytesCnt + TIF.StripByteCounts(ConsolidateCnt);
                if ( BytesCnt >= PlaneBytesCnt ); break; end
            end

            %Consolidate the Strips
            if ( BytesCnt <= PlaneBytesCnt(1) ) && ( ConsolidateCnt > 1 )
                %fprintf(1,'Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
                TIF.StripByteCounts = [BytesCnt; TIF.StripByteCounts(ConsolidateCnt+1:TIF.StripNumber ) ];
                TIF.StripOffsets = TIF.StripOffsets( [1 , ConsolidateCnt+1:TIF.StripNumber] );
                TIF.StripNumber  = 1 + TIF.StripNumber - ConsolidateCnt;
            end
        end
    end

%% read the next IFD address:
    ifd_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);
    %disp(['next ifd at ', num2str(ifd_pos)])

    if isfield( TIF, 'MM_stack' )
        if ( img_last > TIF.MM_stackCnt )
            img_last = TIF.MM_stackCnt;
        end

        %this loop is to read metamorph stacks:
        for ii = img_first:img_last
            TIF.StripCnt = 1;
            %read the image
            fileOffset = PlaneBytesCnt * ( ii - 1 );
            %fileOffset = 0;
            %fileOffset = ftell(TIF.file) - TIF.StripOffsets(1);

            if ( TIF.SamplesPerPixel == 1 )
                IMG.data  = read_plane(fileOffset, IMG.width, IMG.height, 1);
            else
                IMG.red   = read_plane(fileOffset, IMG.width, IMG.height, 1);
                IMG.green = read_plane(fileOffset, IMG.width, IMG.height, 2);
                IMG.blue  = read_plane(fileOffset, IMG.width, IMG.height, 3);
            end

            % print a text timer on the main window, or update the waitbar
            % fprintf(1,'img_read %i img_skip %i\n', img_read, img_skip);
            if exist('waitbar_handle', 'var')
                waitbar( img_read/TIF.MM_stackCnt, waitbar_handle);
            end
            img_read = img_read + 1;
            stack( :,:,ii ) = IMG.data;
        end
        break;
    else
        %this part to read a normal TIFF stack:
        if ( img_skip + 1 >= img_first )
            
            TIF.StripCnt = 1;
            
            %read the image
            if(sum(TIF.StripByteCounts)>2*IMG.width*IMG.height*TIF.BytesPerSample)
                stack = read_full_stack(0, IMG.width, IMG.height,1);   
            else
                if ( TIF.SamplesPerPixel == 1 )
                    IMG.data  = read_plane(0, IMG.width, IMG.height, 1);
                else
                    IMG.red   = read_plane(0, IMG.width, IMG.height, 1);
                    IMG.green = read_plane(0, IMG.width, IMG.height, 2);
                    IMG.blue  = read_plane(0, IMG.width, IMG.height, 3);
                end

                img_read = img_read + 1;
                try
                    stack(:,:, img_read ) = IMG.data;
                catch
                    %stack
                    %IMG
                    error('The file contains dissimilar images: you can only read them one by one');
                end
            end
        else
            img_skip = img_skip + 1;
        end

        if ( img_skip + img_read >= img_last )
            %fprintf(1,'img_skip = %i, img_read = %i, img_last = %i\n', img_skip,img_read,img_last);
            break;
        end
    end
    
end

%% clean-up
fclose(TIF.file);
if exist('waitbar_handle', 'var')
    delete( waitbar_handle );
    clear ('waitbar_handle');
end
drawnow;
%return empty array if nothing was read
if ~ exist( 'stack', 'var')
    stack = [];
end
clear('consolidateStrips','pathname','img_first','img_last','TIF',...
        'img_skip','byte_order','tiff_id','ifd_pos','IMG','num_entries',...
        'file_pos','entry','PlaneBytesCnt','consolidateStrips','BytesCnt',...
        'ConsolidateCnt','i','ii','fileOffset','entry');
clear('IMG');
return;


%============================================================================

function stack = read_full_stack(offset, width, height,planeCnt)
global TIF;

if ( TIF.BitsPerSample(planeCnt) == 0 )
    stack=[];
    return;
end

%determine the type needed to store the pixel values:
switch( TIF.SampleFormat )
    case 1
        classname = sprintf('uint%i', TIF.BitsPerSample(planeCnt));
    case 2
        classname = sprintf('int%i', TIF.BitsPerSample(planeCnt));
    case 3
        if ( TIF.BitsPerSample(planeCnt) == 32 )
            classname = 'single';
        else
            classname = 'double';
        end
    otherwise
        error('unsuported TIFF sample format %i', TIF.SampleFormat);
end

% Preallocate a matrix to hold the sample data:
nz = floor(sum(TIF.StripByteCounts)/(width*height*TIF.BytesPerSample));
stack = zeros(width, height, nz,classname);

% Read the strips and concatenate them:
line = 0; plane = 1;

while ( TIF.StripCnt <= TIF.StripNumber )

    strip = read_strip(offset, width, planeCnt, TIF.StripCnt, classname);
    TIF.StripCnt = TIF.StripCnt + 1;

    % copy the strip onto the data
    if(line+size(strip,2) >= height)
        stack(:, line+1:height,plane) = strip(:,1:height - line);
        stack(:, 1:size(strip,2)-height + line,plane+1) = strip(:,height - line+1:size(strip,2));
        line = size(strip,2)-height + line;
        
        if plane == nz
            if sum(sum(stack(:,:,size(stack,3)))) == 0 
                stack = stack(:,:,1:size(stack,3)-1); 
            end    
            break;
        else
            plane = plane +1;
        end    
    else
        stack(:, line+1:line+size(strip,2),plane) = strip(:,:);
        line = line+size(strip,2);
    end
   
end

% transpose the image (otherwise display is rotated in matlab)
stack = permute(stack,[2,1,3]);

clear('strip','line','plane','classname','nz','TIF');

return;



function plane = read_plane(offset, width, height, planeCnt)

global TIF;

%return an empty array if the sample format has zero bits
if ( TIF.BitsPerSample(planeCnt) == 0 )
    plane=[];
    return;
end

%fprintf(1,'reading plane %i size %i %i\n', planeCnt, width, height);

%determine the type needed to store the pixel values:
switch( TIF.SampleFormat )
    case 1
        classname = sprintf('uint%i', TIF.BitsPerSample(planeCnt));
    case 2
        classname = sprintf('int%i', TIF.BitsPerSample(planeCnt));
    case 3
        if ( TIF.BitsPerSample(planeCnt) == 32 )
            classname = 'single';
        else
            classname = 'double';
        end
    otherwise
        error('unsuported TIFF sample format %i', TIF.SampleFormat);
end

% Preallocate a matrix to hold the sample data:
plane = zeros(width, height, classname);

% Read the strips and concatenate them:
line = 1;
while ( TIF.StripCnt <= TIF.StripNumber )

    strip = read_strip(offset, width, planeCnt, TIF.StripCnt, classname);
    TIF.StripCnt = TIF.StripCnt + 1;

    % copy the strip onto the data
    plane(:, line:(line+size(strip,2)-1)) = strip;

    line = line + size(strip,2);
    if ( line > height )
        break;
    end

end

% Extract valid part of data if needed
if ~all(size(plane) == [width height]),
    plane = plane(1:width, 1:height);
    disp('Cropping data: more bytes read than needed...');
end

% transpose the image (otherwise display is rotated in matlab)
plane = plane';

clear('strip','line','classname','TIF');
return;


%=================== sub-functions to read a strip ===================

function strip = read_strip(offset, width, planeCnt, stripCnt, classname)

global TIF;

%fprintf(1,'reading strip at position %i\n',TIF.StripOffsets(stripCnt) + offset);
StripLength = TIF.StripByteCounts(stripCnt) ./ TIF.BytesPerSample(planeCnt);

%fprintf(1, 'reading strip %i\n', stripCnt);
fseek(TIF.file, TIF.StripOffsets(stripCnt) + offset, 'bof');
bytes = fread( TIF.file, StripLength, classname, TIF.BOS );

if ( length(bytes) ~= StripLength )
    error('End of file reached unexpectedly.');
end

strip = reshape(bytes, width, StripLength / width);
clear('bytes','StripLength','TIF');
return;


%===================sub-functions that reads an IFD entry:===================


function [nbBytes, matlabType] = convertType(tiffType)
switch (tiffType)
    case 1
        nbBytes=1;
        matlabType='uint8';
    case 2
        nbBytes=1;
        matlabType='uchar';
    case 3
        nbBytes=2;
        matlabType='uint16';
    case 4
        nbBytes=4;
        matlabType='uint32';
    case 5
        nbBytes=8;
        matlabType='uint32';
    case 11
        nbBytes=4;
        matlabType='float32';
    case 12
        nbBytes=8;
        matlabType='float64';
    otherwise
        error('tiff type %i not supported', tiffType)
end
return;

%===================sub-functions that reads an IFD entry:===================

function  entry = readIFDentry()

global TIF;
entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.BOS);
entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.BOS);
%disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);

[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32', TIF.BOS);
    %disp(strcat('offset = ', num2str(offset)));
    fseek(TIF.file, offset, -1);
end

if TIF.entry_tag == 33629   %special metamorph 'rationals'
    entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.BOS);
else
    if entry.tiffType == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);
    else
        entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.BOS);
    end
end
if ( entry.tiffType == 2 ); entry.val = char(entry.val'); end
clear('offset','TIF');
return;


%==============distribute the metamorph infos to each frame:
function [info, stack, wavelength, private2 ] = extractMetamorphData(imgCnt)

global TIF;

info = [];
stack = [];
wavelength = [];
private2 = [];

if TIF.MM_stackCnt == 1
    return;
end

left  = imgCnt - 1;

if isfield( TIF, 'info' )
    S = length(TIF.info) / TIF.MM_stackCnt;
    info = TIF.info(round(S*left+1):round(S*left+S));
end

if isfield( TIF, 'MM_stack' )
    S = length(TIF.MM_stack) / TIF.MM_stackCnt;
    stack = TIF.MM_stack(round(S*left+1):round(S*left+S));
end

if isfield( TIF, 'MM_wavelength' )
    S = length(TIF.MM_wavelength) / TIF.MM_stackCnt;
    wavelength = TIF.MM_wavelength(round(S*left+1):round(S*left+S));
end

if isfield( TIF, 'MM_private2' )
    S = length(TIF.MM_private2) / TIF.MM_stackCnt;
    private2 = TIF.MM_private2(S*left+1:S*left+S);
end

clear('left','S','TIF');
return;
