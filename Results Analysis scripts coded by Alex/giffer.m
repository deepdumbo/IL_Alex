%{
    last updated on 2020/03/06, Hossein

    Given a stack of pictures, generate a gif. 
%}

%load('C:\Users\Alex Ensworth\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Nov 8 2019\Saved_MAT_files\Axial_Mask_erosion_1_MEDI RG_PDF_TKD.mat')
%load('C:\Users\Alex Ensworth\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Nov 8 2019\Saved_preQSM_data\20191108mrdev3126917723.mat')

%PhaseData=double(dataScan.phaseReconstructed);
QSMdat=images.PrettyQSM;
matrix_size=size(QSMdat);

%Filename that the gif will be saved as:
filename = 'QSMTKD.gif';

fig=figure;
for i=1:matrix_size(3)
    imagesc(QSMdat(:,:,i),[-0.3,0.3]);
    axis image;
    colorbar
    colormap(bone)
    drawnow
    frame=getframe(fig);
    im=frame2im(frame);

    [A,map]=rgb2ind(im,256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.035);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.035);
    end
end