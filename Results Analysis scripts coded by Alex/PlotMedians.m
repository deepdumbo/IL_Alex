FigSaveLocation='D:\AlexEnsworth\CurrentDataAndFigs\Figures\';

MidVal=zeros(1,size(images.iFreq,3));
for i=1:size(images.iFreq,3)
    MidVal(i)=median(nonzeros(images.iFreq(:,:,i)));
end

MidVal(isnan(MidVal))=0;

plot(MidVal)

saveas(gcf, [FigSaveLocation, images.dataset '_' images.UnwrapType '_' images.TU 'combined_echoes'])
print(gcf, [FigSaveLocation, images.dataset '_' images.UnwrapType '_' images.TU 'combined_echoes'],'-dpng','-r600')