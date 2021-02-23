%Mehdi Ibrahimli
%ID 2040467
clear;
clc;
%read the hdf
data=hdfread('myfile.L1R','EO1H1660322004097110PX.L1R');

%random band and area
clipdata=data(1:256,50,:);

%data 8 bit transformation and histogram stretching and equalization
clipdata=(clipdata/256);
m=max(max(clipdata));
clipdata=(510/m)*clipdata-255;
clipdata=uint8(clipdata);
clipdata=histeq(clipdata);

% show image
clipdata=permute(clipdata,[1,3,2]);
colormap(gray);
image(clipdata);
%% DFT
FFT = fft2(clipdata);
FD=fftshift(log(1+(abs(FFT))));
imagesc(FD);

%high pass filter
FDH = fftshift(FFT);
FDH(104:154,104:154)=0;
FDH = ifftshift(FDH);
FDH = ifft2(FDH);
imagesc(FDH,[-150,150]);

%low pass filter
FDL=fftshift(FFT);
FDL=FDL-FDH;
FDL = ifftshift(FDL);
FDL = abs(ifft2(FDL));
imagesc(FDL,[0,250]);





