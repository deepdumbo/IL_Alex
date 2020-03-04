clear all;clc;
addpath('LBV');
addpath('PDF');
addpath('SHARP');
addpath('RESHARP');
addpath('../Common');
load ../human_RDF;
load Mask_BGRM
Mask = Mask_BGRM;
radius = 4;
M1 = SMV(Mask, matrix_size, voxel_size, radius)>0.999;
RDF_LBV = (LBV(iFreq, Mask, matrix_size, voxel_size, 1e-4,4,1)+0.05).*M1;
RDF_SHARP = SHARP(iFreq, Mask, matrix_size, voxel_size, radius,0.03).*M1;
RDF_RESHARP = RESHARP(iFreq, Mask, matrix_size, voxel_size, radius, 0.01).*M1;
RDF_PDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir).*M1;
z = 69;
gray_range = [-0.2 0.2];
close all;
figure;imshow(RDF_LBV(:,:,z)',gray_range);
figure;imshow(RDF_SHARP(:,:,z)',gray_range);
figure;imshow(RDF_RESHARP(:,:,z)',gray_range);
figure;imshow(RDF_PDF(:,:,z)',gray_range);
