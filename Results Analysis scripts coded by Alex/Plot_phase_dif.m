% clear
% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end

set(gcf, 'Position', get(0,'Screensize'))
% for i=1:size(PhaseData,4)
%     subplot(2,3,i)
%     imagesc(PhaseData(:,:,65,i))%,[0,15]);
%     axis image
%     title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, i));%, PDF_TKD.RDFtype, images.slice))
%     colorbar
%     impixelinfo
% end

PhaseData=images.PhaseData;
iFreq=images.iFreq;
    subplot(2,3,1)
%     imagesc(PhaseData(:,:,45,1),[-7,3]);
    imagesc(PhaseData(:,:,images.slice,1));
%     imagesc(iFreq_raw(:,:,45,1))%,[-4,4]);
    axis image
    title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, 1));%, PDF_TKD.RDFtype, images.slice))
% title(sprintf('Wrap phase Echo number %d', 1));
    colorbar
    impixelinfo
    
    subplot(2,3,2)
%     imagesc(PhaseData(:,:,45,2),[-9,9]);
    imagesc(PhaseData(:,:,images.slice,2));
%     imagesc(iFreq_raw(:,:,45,2))%,[-4,4]);    
    axis image
    title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, 2));%, PDF_TKD.RDFtype, images.slice))
%     title(sprintf('Wrap phase Echo number %d', 2));
    colorbar
    impixelinfo
    
    subplot(2,3,3)
%     imagesc(PhaseData(:,:,45,3),[-10,13]);
    imagesc(PhaseData(:,:,images.slice,3));
%     imagesc(iFreq_raw(:,:,45,3))%,[-4,4]);    
    axis image
%     title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, 3));%, PDF_TKD.RDFtype, images.slice))
    title(sprintf('Wrap phase Echo number %d', 3));
    colorbar
    impixelinfo
    
    subplot(2,3,4)
%     imagesc(PhaseData(:,:,45,4),[-10,18]);
    imagesc(PhaseData(:,:,images.slice,4));
%     imagesc(iFreq_raw(:,:,45,4))%,[-4,4]);    
    axis image
%     title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, 4));%, PDF_TKD.RDFtype, images.slice))
    title(sprintf('Wrap phase Echo number %d', 4));
    colorbar
    impixelinfo
    
    subplot(2,3,5)
%     imagesc(PhaseData(:,:,45,5),[-15,23]);
    imagesc(PhaseData(:,:,images.slice,5));
%     imagesc(iFreq_raw(:,:,45,5))%,[-4,4]);    
    axis image
    title(sprintf('Phase Unwrap using %s, Echo number %d', images.UnwrapType, 5));%, PDF_TKD.RDFtype, images.slice))
%     title(sprintf('Wrap phase Echo number %d', 5));
    colorbar
    impixelinfo

    
subplot(2,3,6)
imagesc(iFreq(:,:,images.slice));%,[-2,6]);
axis image
title(sprintf('The combined, unwrapped phase data; %s', images.UnwrapType));%, PDF_TKD.RDFtype, images.slice))
colorbar
impixelinfo


% for i=1:size(iFreq_raw,4)
%     subplot(2,3,i)
%     imagesc(iFreq_raw(:,:,65,i));%,[-0.3,0.3]);
%     axis image
%     title(sprintf('Phase Wrapped data, Echo number %d', i));%, PDF_TKD.RDFtype, images.slice))
%     colorbar
%     impixelinfo
% end



% colormap(bone)

colormap(gray)

saveas(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_echo_phase_unwrap'])
print(gcf, ['./Plots/', images.dataset '_' images.UnwrapType '_echo_phase_unwrap'],'-dpng','-r600')

% saveas(gcf, ['./Plots/', images.dataset '_echo_phase_wrapped'])
% print(gcf, ['./Plots/', images.dataset '_echo_phase_wrapped'],'-dpng','-r600')
