function If = DOGfilter(I,s1,s2,padsize,kernsize)
%Performs a difference of gaussian filter on an image or z-stack I
%s1: small gaussian sigma (suggested for single spots: 1)
%s2: large gaussian sigma (suggested for single spots: ~twice the PSF sigma)

%optional arguments (enter empty array to set to default)
%padsize: size of padding (in pixels) added to the image; default = 30;
%kernsize: size of convolution kernel; default = 21;

if isempty(padsize)
    padsize = 30;
end

if isempty(kernsize)
    kernsize = 21;
end

g1 = fspecial('Gaussian', kernsize, s1);
g2 = fspecial('Gaussian', kernsize, s2);
dog = g1 - g2;

If = zeros(size(I));

if ndims(I) ==3
    [nx,ny,nz] = size(I);    
    for i=1:nz
        B = padarray(I(:,:,i),[padsize,padsize],'symmetric');
        Idog = conv2(double(B), dog, 'same');
        If(:,:,i) = Idog(padsize+1:padsize+nx,padsize+1:padsize+ny);    
    end    
else
    [nx,ny] = size(I);
    B = padarray(I,[padsize,padsize],'symmetric');
    Idog = conv2(double(B), dog, 'same');
    If = Idog(padsize+1:padsize+nx,padsize+1:padsize+ny);
end
end

