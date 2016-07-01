clc
clear all
bb=8; % block size
K=256; % number of atoms in the dictionary
img =imread('Test_Fig2_Missing.png');

[N,M,dim]=size(img);
img = double(img);
%Compute mask and extracting its patches
Mask = double(~(img(:,:,1)==0));
blkMask=im2col(Mask,[bb,bb],'sliding');  % distinct  sliding
img_yuv = rgb2ycbcr(uint8(img));
img_inpaint_yuv = zeros(size(img_yuv));
% Interpolation CbCr Componet
img_inpaint_yuv(:,:,2) = Interpolation(double(img_yuv(:,:,2)),~Mask); 
img_inpaint_yuv(:,:,3) = Interpolation(double(img_yuv(:,:,3)),~Mask);
IMin0 = double(img_yuv(:,:,1)); 

load Dict
load Coeff

% Creating the output image
imag_Y=ImageRecover(IMin0,Dict,Coeff); 
imag_Y=max(min(imag_Y,255),0);

img_inpaint_yuv(:,:,1) = imag_Y; 
img_inpaint_rgb = ycbcr2rgb(uint8(img_inpaint_yuv));
imshow(uint8(img_inpaint_rgb))
imwrite(img_inpaint_rgb,strcat('KSVD_Result_','iter_25','.png'),'png')