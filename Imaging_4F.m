%% 4F imaging system using fourier optics


clear all;
close all;
clc;


a=SequentialOpticalModel;
%% set the main parameters of all matrices and filters used

dim=512;      % the dimension of the array image or filter
center= dim/2 +1 ;
s= 3;
f=1/(dim*s); % frequency spacing 

%% reading the input image

img = imread('InputImage.jpg');
img = imresize(img,[dim,dim]);                                 
                                      
                                 
                                      
                                      
%% Fourier transform of the Image
FT_img= fftshift(fft2(fftshift(img(:,:,3))));
ftimg_H_SingleSlit = FT_img.*a.H_SingleSlit(dim,center);
ftimg_horizontal_double_slit = FT_img.*a.horizontal_double_slit(dim,center);

ftimg_vertical_single_slit = FT_img.*a.vertical_single_slit(dim,center);
ftimg_vertical_double_slit = FT_img.*a.vertical_double_slit(dim,center);

ft_pinhole = FT_img.*a.pinhole_filter(dim,f);



% calculate inverse FT and the absolute value to plot it  
img1_H_SingleSlit = abs(fftshift(ifft2(fftshift(ftimg_H_SingleSlit))));
img2_horizontal_double_slit = abs(fftshift(ifft2(fftshift(ftimg_horizontal_double_slit))));

img3_vertical_single_slit = abs(fftshift(ifft2(fftshift(ftimg_vertical_single_slit))));
img4_vertical_double_slit = abs(fftshift(ifft2(fftshift(ftimg_vertical_double_slit))));
img_pinhole = abs(fftshift(ifft2(fftshift(ft_pinhole))));



%%    >>>>>> PLOTTING



a.draw("pinhole",a.pinhole_filter(dim,center),img_pinhole)
a.draw('Horizontal Single slit',a.H_SingleSlit(dim,center),img1_H_SingleSlit)
a.draw('Horizontal double slit',a.horizontal_double_slit(dim,center),img2_horizontal_double_slit)
a.draw('verical single slit',a.vertical_single_slit(dim,center),img3_vertical_single_slit)
a.draw('Horizontal Single slit',a.vertical_double_slit(dim,center),img4_vertical_double_slit)








