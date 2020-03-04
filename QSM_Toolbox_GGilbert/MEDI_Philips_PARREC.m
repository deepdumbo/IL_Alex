function MEDI_Philips_PARREC()
filesPath=strcat(pwd,'/QSM_toolbox');
addpath(genpath(filesPath));
%**************************************************************************
% Path to input PAR/REC files (magnitude AND phase)
%**************************************************************************
[FileName,PathName] = uigetfile({'*.*',  'All Files (*.*)'},'Please select the PAR/REC file');
FileName = strtok(FileName, '.');

%**************************************************************************
% We read the input data
%**************************************************************************
% In the case of R2/R3/R4 data
[Data info]=ParRecReader(PathName,FileName,1);
% In the case of R5 data
%[Data info]=ParRecReaderR5(PathName,FileName,1);
[Nx,Ny,NSlices,tt] = size(Data);


%**************************************************************************
% In the case of highly anisotropic data, we may want to perform some
% interpolation (also to speed-up the calculation)
%**************************************************************************
interpolation=0.5;
for w=1:info.NEchos;
    for k=1:NSlices
        Data2(:,:,k,w)=imresize(Data(:,:,k,w),interpolation);
    end
end
Data=Data2;
clear Data2;

info.xSpace=info.xSpace/interpolation;
info.ySpace=info.ySpace/interpolation;

%**************************************************************************
% We can also add some zero padding if we want
%**************************************************************************
padding_pixels=[0 0 8 0];

Data=padarray(Data,padding_pixels,'post');

[Nx,Ny,NSlices,tt] = size(Data);
MagnData=abs(Data);
PhaseData=angle(Data);
%**************************************************************************
% Loop over the different echoes
%**************************************************************************
for jj=1:info.NEchos;
    
    %**********************************************************************
    % Laplacian unwrapping (Li et al.)
    %**********************************************************************
    display('Phase unwrapping');
    [UPhaseData]=LaplacianUnwrap(PhaseData(:,:,:,jj),[info.ySpace info.xSpace],(info.zThick+info.zSpace));
    PhaseData(:,:,:,jj)=UPhaseData; 
        
 end
clear UPhaseData;

%*************************************************************************
%Combination of the phase information using a weighted non-linear regression (Gilbert
%et al.)
% *************************************************************************
if (info.NEchos>1)
    display('Combination of the echoes');
    [iFreq N_std]=CombinationEcho(MagnData,PhaseData,info);
    delta_TE=(info.TE(2)-info.TE(1));
else
    iFreq=PhaseData./info.TE; % Normalisation
    N_std=1./MagnData;
    delta_TE=info.TE;
end
N_std(isnan(N_std))=0;
N_std(isinf(N_std))=0;

%**************************************************************************
% Calculation of a (rough) brain mask. Which echo to use?
%**************************************************************************
Mask=logical(AutoMaskRec(MagnData(:,:,:,2))); 

%**************************************************************************
% Calculation of the local field using the RESHARP method
%**************************************************************************
[RDF Mask]=RESHARP(iFreq,Mask,size(iFreq),[info.xSpace info.ySpace (info.zSpace+info.zThick)],5,0.01);
N_std=N_std.*Mask;
RDF=RDF*delta_TE;

% %**************************************************************************
% % We save the RDF image 
% %**************************************************************************
Image_Write=RDF; 
for k=1:size(Image_Write,3)
    Image_Write(:,:,k)=fliplr(Image_Write(:,:,k));
end
% hdr = analyzeWrite(Image_Write,[PathName,'\',FileName,'RDF'], [info.ySpace info.xSpace (info.zThick+info.zSpace)], [], 0);
hdr = analyzeWrite(Image_Write,[PathName,FileName,'RDF'], [info.ySpace info.xSpace (info.zThick+info.zSpace)], [], 0);


%**************************************************************************
% We save the data, to be used by the MEDI_linear code
%**************************************************************************
%We evaluate the B0 direction (NOTE: we assume an axial scan and a single stack!!!)
% iop=info{1}.ImageOrientationPatient;
% r=iop(1:3);  c=iop(4:6); s=cross(r',c');
% R = [r(1) c(1) s(1); r(2) c(2) s(2); r(3) c(3) s(3);];
Rx=[1 0 0; 0 cos(info.rotX) -sin(info.rotX); 0 sin(info.rotX) cos(info.rotX)];
Ry=[cos(info.rotY) 0 sin(info.rotY); 0 1 0; -sin(info.rotY) 0 cos(info.rotY)];
Rz=[cos(info.rotZ) -sin(info.rotZ) 0; sin(info.rotZ) cos(info.rotZ) 0; 0 0 1];
Sys=Rx*Ry*Rz; clear Rx Ry Rz;
B0_dir=Sys*([0 0 1].'); 
CF=127.72e6;
iMag=sum(sqrt(MagnData.^2),4).*Mask; % Which echo should we use?
matrix_size=size(iFreq);
voxel_size=[info.xSpace info.ySpace (info.zSpace+info.zThick)];

save('Philips.mat','B0_dir','CF','Mask','iFreq','iMag','delta_TE','matrix_size','voxel_size','RDF','N_std');

%**************************************************************************
% We perform the dipole inversion
%**************************************************************************
QSM=MEDI_linear('lambda',600,'merit',1).*Mask;

%**************************************************************************
% We save the QSM image 
%**************************************************************************
Image_Write=QSM;
for k=1:size(Image_Write,3)
    Image_Write(:,:,k)=fliplr(Image_Write(:,:,k));
end
% hdr = analyzeWrite(Image_Write,[PathName,'\',FileName,'MEDI'], [info.ySpace info.xSpace (info.zThick+info.zSpace)], [], 0);
hdr = analyzeWrite(Image_Write,[PathName,FileName,'MEDI'], [info.ySpace info.xSpace (info.zThick+info.zSpace)], [], 0);


