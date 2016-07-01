% Author: Chongyang
% Date: 2016/06/05
% Main function for image inpainting based on OMP and KSVD
% Case 1: Interpolation
% Case 2: OMP
% Case 3: KSVD
clc
clear all;
close all;


bb=8; % block size
K=256; % number of atoms in the dictionary
overlap = 3; %number of pixels between consecutive patches
J = 10;  % Number of iteration in KSVD 

img =imread('Test_Fig2_Missing.png');      

Case = 2;
switch Case
    %% Interpolation
    case 1 
        Mask = double(~(img(:,:,1)==0));
        for channel = 1:3
            IMin0 = img(:,:,channel);
            IMin0= double(IMin0);
            out(:,:,channel) = Interpolation(IMin0,~Mask);
        end
        imshow(uint8(out))
        imwrite(uint8(out),strcat('Test1_Interpolation_Result','.png'),'png')
    %% OMP
    case 2 
        [N,M,dim]=size(img);
        img = double(img);
        NN = ceil((N-bb)/overlap) * overlap + bb;
        MM = ceil((M-bb)/overlap) * overlap + bb;
        img_new = 255 * zeros(NN,MM,dim);
        img_new(1:N,1:M,:) = img;  
        %imshow(uint8(img_new))

        
        %Compute mask and extracting its patches
        Mask = double(~(img_new(:,:,1)==0));
        blkMask = overlap_im2col(Mask, bb, overlap);
        
        
        img_yuv = rgb2ycbcr(uint8(img_new));
        img_inpaint_yuv = zeros(size(img_yuv));
        % Interpolation CbCr Componet
        img_inpaint_yuv(:,:,2) = Interpolation(double(img_yuv(:,:,2)),~Mask); 
        img_inpaint_yuv(:,:,3) = Interpolation(double(img_yuv(:,:,3)),~Mask);


        IMin0 = double(img_yuv(:,:,1)); 

        % The dictionary DCT
        DCT=zeros(bb,sqrt(K));
        for k=0:1:sqrt(K)-1,
            V=cos([0:1:bb-1]'*k*pi/sqrt(K));
            if k>0 
                V=V-mean(V); 
            end;
            DCT(:,k+1)=V/norm(V);
        end;
        DCT=kron(DCT,DCT);

        % Extracting the noisy image patches
        blkMatrixIm = overlap_im2col(IMin0, bb, overlap);
        
        
     
        sigma = 0.001;  % 0.005  0.01
        rc_min = 0.01; % rc_min: minimal residual correlation before stopping pursuit
        max_coeff = 10; % max_coeff: sparsity constraint for signal representation 10

        % Inpainting the Patches (OMP)
        Coeff=OMP_Inpainting(DCT,blkMatrixIm.*blkMask,blkMask,sigma,rc_min,max_coeff); %256*146405
        Coeff = full(Coeff);
        
        
        % Creating the output image
        imag_Y = overlap_col2im(DCT*Coeff, blkMask, bb, overlap, size(img_new));
        imag_Y=max(min(imag_Y,255),0);

        img_inpaint_yuv(:,:,1) = imag_Y; 
        img_inpaint_rgb = ycbcr2rgb(uint8(img_inpaint_yuv));
        img_out = img_inpaint_rgb(1:N,1:M,:);
        imshow(uint8(img_out))
        imwrite(img_out,strcat('Test2_OMP_Result_',num2str(sigma),'_.png'),'png')
   %% KSVD
    case 3
        [N,M,dim]=size(img);
        img = double(img);
        NN = ceil((N-bb)/overlap) * overlap + bb;
        MM = ceil((M-bb)/overlap) * overlap + bb;
        img_new = 255 * zeros(NN,MM,dim);
        img_new(1:N,1:M,:) = img;  
        %imshow(uint8(img_new))
    
        
        %Compute mask and extracting its patches
        Mask = double(~(img_new(:,:,1)==0));
        blkMask = overlap_im2col(Mask, bb, overlap);
        
        
        
        img_yuv = rgb2ycbcr(uint8(img_new));
        img_inpaint_yuv = zeros(size(img_yuv));
        % CbCr Componet Interpolation
        img_inpaint_yuv(:,:,2) = Interpolation(double(img_yuv(:,:,2)),~Mask); 
        img_inpaint_yuv(:,:,3) = Interpolation(double(img_yuv(:,:,3)),~Mask);


        IMin0 = double(img_yuv(:,:,1)); 

        % The dictionary DCT
        DCT=zeros(bb,sqrt(K));
        for k=0:1:sqrt(K)-1,
            V=cos([0:1:bb-1]'*k*pi/sqrt(K));
            if k>0 
                V=V-mean(V); 
            end;
            DCT(:,k+1)=V/norm(V);
        end;
        DCT=kron(DCT,DCT);

        % Extracting the noisy image patches
        blkMatrixIm = overlap_im2col(IMin0, bb, overlap);
        
        sigma = 0.01;  % 0.005 0.01
        rc_min = 0.01; % rc_min: minimal residual correlation before stopping pursuit
        max_coeff = 10; % max_coeff: sparsity constraint for signal representation 10

        
        % Inpainting the Patches (K-SVD)
        [Dict Coeff]=KSVD_Inpainting(DCT,blkMatrixIm,blkMask,sigma,rc_min,max_coeff,J);  %DCT

      
        % Creating the output image
        imag_Y = overlap_col2im(Dict*Coeff, blkMask, bb, overlap, size(img_new));
        imag_Y=max(min(imag_Y,255),0);

        img_inpaint_yuv(:,:,1) = imag_Y; 
        img_inpaint_rgb = ycbcr2rgb(uint8(img_inpaint_yuv));
        img_out = img_inpaint_rgb(1:N,1:M,:);
        imshow(uint8(img_out))
        imwrite(img_out,strcat('Test3_KSVD_Result_',num2str(sigma),'_.png'),'png')
        
end
