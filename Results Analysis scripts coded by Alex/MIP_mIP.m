%{
   last updated on 2020/03/06   

   when given a stack of slices from a QSM image, this script merges the stack to generate one slice. 
   the final slice is made by selecting only one voxel from the coloumn of voxels in the stack. 
   the voxel selection is based on either min or max signal. 
   
   input: QSM images, slice range
   output: figure showing average QSM slice, the slice from min signal intensity, slice from max signal intensity
%}

%Maximum and Minimum intensity projections (MIPs and mIPs)

%% Load the file

[file,path] = uigetfile('*.mat');
if isequal(file,0)
   disp('User selected Cancel');
   return
else
   load(fullfile(path,file));
end

%% User defined slice selection

slice_start = 72;
slice_stop = 86;

mid_slice=ceil((slice_start + slice_stop)/2);

%% Get the mIP and MIP

subset_QSM=images.QSM(:,:,slice_start:slice_stop); %this step makes the max and min calc easy

MIP = max(subset_QSM,[],3); %Take the maximum value of each pixel, through the 3rd dimension
mIP = min(subset_QSM,[],3); %Same, but the minimum value

%% Plot

set(gcf, 'Position', get(0,'Screensize'))

% The midslice QSM image
subplot (1,3,1)
imagesc(images.QSM(:,:,mid_slice),[-0.3,0.3])
axis image
title(sprintf('%s QSM, slice %d', images.QSMtype, mid_slice))
colorbar
impixelinfo

% The maximum intensity projection
subplot (1,3,2)
imagesc(MIP(:,:),[-0.1,0.4])
axis image
title(sprintf('MIP, slices %d to %d', slice_start, slice_stop))
colorbar
impixelinfo

% The minimum intensity projection
subplot (1,3,3)
imagesc(mIP(:,:),[-0.4,0.1])
axis image
title(sprintf('mIP, slices %d to %d', slice_start, slice_stop))
colorbar
impixelinfo

colormap(bone)

%% Save the plots
saveas(gcf, ['./MIP_mIP_Plots/' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_MIP_mIP_Slices_' num2str(slice_start) 'to' num2str(slice_stop)])
print(gcf, ['./MIP_mIP_Plots/' images.UnwrapType '_' images.RDFtype '_' images.QSMtype '_MIP_mIP_Slices_' num2str(slice_start) 'to' num2str(slice_stop)],'-dpng','-r600')
