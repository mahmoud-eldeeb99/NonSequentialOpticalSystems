%% 4F imaging system using fourier optics


clear all;
close all;
clc;


a=NonSequentialOpticalModel;
%% set the main parameters of all matrices and filters used

dim=512;      % the dimension of the array image or filter
center= dim/2 +1 ;
s= 3;
f=1/(dim*s); % frequency spacing 

%% reading the input image

img = imread('InputImage.jpg');
img = imresize(img,[dim,dim]);

%bin_object = rgb2gray(object);

%% creatig the filters

% 1.Horizontal Single Slit 
H_SingleSlit= a.H_SingleSlit(dim,center);

% 2.horizontal double slit 
horizontal_double_slit=a.horizontal_double_slit(dim,center);
              

%% 3.Vertical Single Slit 
vertical_single_slit=a.vertical_single_slit(dim,center);
                                      
                                      
% 4.Vertical Double Slit
vertical_double_slit=a.vertical_double_slit(dim,center);
                                      
                                      
%% 5.pinhole                                       
                                      
pinhole_filter=a.pinhole_filter(dim,f);                                  
                                      
                                      
%% Fourier transform of the Image
FT_img= fftshift(fft2(fftshift(img(:,:,3))));
ftimg_H_SingleSlit = FT_img.*H_SingleSlit;
ftimg_horizontal_double_slit = FT_img.*horizontal_double_slit;

ftimg_vertical_single_slit = FT_img.*vertical_single_slit;
ftimg_vertical_double_slit = FT_img.*vertical_double_slit;

ft_pinhole = FT_img.*pinhole_filter;



% calculate inverse FT and the absolute value to plot it  
img1_H_SingleSlit = abs(fftshift(ifft2(fftshift(ftimg_H_SingleSlit))));
img2_horizontal_double_slit = abs(fftshift(ifft2(fftshift(ftimg_horizontal_double_slit))));

img3_vertical_single_slit = abs(fftshift(ifft2(fftshift(ftimg_vertical_single_slit))));
img4_vertical_double_slit = abs(fftshift(ifft2(fftshift(ftimg_vertical_double_slit))));
img_pinhole = abs(fftshift(ifft2(fftshift(ft_pinhole))));



%%    >>>>>> PLOTTING



%% 1. horezontal single slit 

figure('Name', 'Horizontal Single slit');
%set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('summer');

subplot(1,3,1);
imagesc(H_SingleSlit);axis('image');
title('Horozental single slit');

subplot(1,3,2);
imagesc(img1_H_SingleSlit);axis('image');
title('horzental single sit');

subplot(2,3,3);
mesh(img1_H_SingleSlit); 
title('Intensity profile');



%% 2.horizontal double slit 
figure('Name',  'Horizontal double slit');
colormap('summer');
subplot(1,3,1);
imagesc(horizontal_double_slit);axis('image');
title('Horizontal double slit');

subplot(1,3,2);
imagesc(img2_horizontal_double_slit);axis('image');
title('horizontal double slit image');


subplot(2,3,3);
mesh(img2_horizontal_double_slit); 
title('Intensity profile');


%% 3. verical single slit 

figure('Name', 'verical single slit');
colormap('summer');

subplot(1,3,1);
imagesc(vertical_single_slit);axis('image');
title('verical single slit filter');

subplot(1,3,2);
imagesc(img3_vertical_single_slit);axis('image');
title('vertical single sit image');

subplot(2,3,3);
mesh(img3_vertical_single_slit); 
title('Intensity profile');

%% 4. vertical double slit 

figure('Name', 'Horizontal Single slit');
%set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('summer');

subplot(1,3,1);
imagesc(vertical_double_slit);axis('image');
title('vertical double slit filter');

subplot(1,3,2);
imagesc(img4_vertical_double_slit);axis('image');
title('vertical double slit image');

subplot(2,3,3);
mesh(img4_vertical_double_slit); 
title('Intensity profile');



%% 5.pinhole
figure('Name', 'pinhole');
colormap('summer');
subplot(1,3,1);
imagesc(pinhole_filter);axis('image');
title('pinhole filter');

subplot(1,3,2);
imagesc(img_pinhole);axis('image');
title('pinhole Image');

subplot(2,3,3);
mesh(img_pinhole); 
title('Intensity profile');

                 
