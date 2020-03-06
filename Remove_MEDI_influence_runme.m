%{
    QSM Analysis main script
    The impact of unwrapping and background removal on QSM is explored. 
    
    Input:
        Pre-QSM mat files containing magnitue and phase maps
    Outputs:
        mat file with specific MaskErode, UnwrapType, RDFtype, QSMtype

    The sections of this script perform the following tasks:

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time=clock;
images.start_time=start_time;

%% Data download

% *************************************************************************
% Make sure the current directory is 'Ives Levesque, Mr - mtbi_qmri'
% *************************************************************************
cd 'C:\Users\alexe\McGill University\Ives Levesque, Mr - mtbi_qmri'

% *************************************************************************
% Save location on my hard drive
% *************************************************************************
DataSaveLocation='D:\AlexEnsworth\CurrentDataAndFigs\Data\';
FigSaveLocation='D:\AlexEnsworth\CurrentDataAndFigs\Figures\';

% *************************************************************************
% Save location on Ives' OneDrive:
% *************************************************************************
% DataSaveLocation='./data/Analysed/Saved_MAT_files/';
% FigSaveLocation='./figures/Plots/';

%%
% *************************************************************************
% Select the appropriate set of data
%{
    Load mat files that are previosly generated from DICOM files. 
    The matfiles contain 
        dataScan, 
        ScannerType 
    I   Images
    
%}
% *************************************************************************
% Jan 22 data (Siemens):
load('./data/Processed/Saved_preQSM_data/20200122mrdev158941777scale_phase.mat')
ScannerType='Siemens';
images.dataset='Siemens_Jan22';

% Dec 13 Run 1:************************************************************
% load('./data/Processed/Saved_preQSM_data/20191213mrdev3126917723run1.mat')
% ScannerType='Philips';
% images.dataset='Dec_13_run1';

% Dec 13 Run 2:************************************************************
% load('./data/Processed/Saved_preQSM_data/20191213mrdev3126917723run2.mat')
% ScannerType='Philips';
% images.dataset='Dec_13_run2';

 %Nov 8 data:**************************************************************
% load('./data/Processed/Saved_preQSM_data/20191108mrdev3126917723.mat')
% ScannerType='Philips';
% images.dataset='Nov_8';


%% Data preperation
%{
    Data from dataScan is extracted
    Depending on the scanner, scanner-dependent variables are set
    dataScan has the following variables
        "Magnitude", "Phase", "voxelSize", "matrixSize", "CF", "TE" , "B0_dir"
%}
gamma=267.513E6;

if strcmp(ScannerType,'Siemens')
    mask_red_num=104;                  %This is to only have the mask cover the brain, and not the neck. 134 seems to be a good number.
    images.slice=45;
    MagnData=double(dataScan.magnitude);
    PhaseData=double(dataScan.phase);
else
    mask_red_num=154;
    images.slice=65;
    MagnData=double(dataScan.magn_rescaleFP);
    PhaseData=double(dataScan.phaseReconstructed);
end

voxel_size=dataScan.voxelSize;
images.voxel_size=voxel_size;
matrixSize=dataScan.matrixSize;
CF=dataScan.CF;
    images.CF=CF;
mainField=2*pi*CF/gamma; 
    images.mainField=mainField;
TE=dataScan.TE;
    images.TE=TE;
if length(TE)>1
    delta_TE=TE(2)-TE(1);
end
    images.delta_TE=delta_TE;
B0_dir=dataScan.B0_dir';
    images.B0_dir=B0_dir;

iField_dat = MagnData.*exp(1i*PhaseData);


%Axial Orientation - The default

images.orient='Axial';
% Find the odd spot and then reduce it to make it even. Note.. this only works if there is only one odd number
if mod(size(iField_dat,1),2)==0 
    if mod(size(iField_dat,2),2)==0
        if mod(size(iField_dat,3),2)==0
            iField=iField_dat(:,:,:,:);
            PhaseData=PhaseData(:,:,:,:);
            MagnData=MagnData(:,:,:,:);
        else
            iField=iField_dat(:,:,1:(end-1),:);
            PhaseData=PhaseData(:,:,1:(end-1),:);
            MagnData=MagnData(:,:,1:(end-1),:);
        end
    else
        iField=iField_dat(:,1:(end-1),:,:);
        PhaseData=PhaseData(:,1:(end-1),:,:);
        MagnData=MagnData(:,1:(end-1),:,:);
    end
        
else
    iField=iField_dat(1:(end-1),:,:,:);
    PhaseData=PhaseData(1:(end-1),:,:,:);
    MagnData=MagnData(1:(end-1),:,:,:); %The TKD implementation complained that there was an odd number of slices, so we made it even
end

Big_Size=size(iField);
matrix_size=Big_Size(1:3); 
    images.matrix_size=matrix_size;


images.mask_red_num=mask_red_num;
matrix_size_mask=matrix_size;       %This is used along with Mask_red (reduced), and iMag_mask
matrix_size_mask(3)=mask_red_num;                               


images.iField=iField;

num_echo=size(TE,2);

%% Magnitude Image
% Compute magnitude image
iMag = sqrt(sum(abs(iField).^2,4));
iMag_mask=iMag(:,:,1:mask_red_num);
images.iMag=iMag;

%% Mask
% Use FSL BET to extract brain mask
%{
	Raw Imaginary frequency "iFreq_raw" is generated by applying the mask to phase data
	The mask is also applied to iField and is subtracted from it. 

    Error with MEDI_toolbox
        bet.mexw64 not found eventhough this files is added to path
    Hossein, 2020/03/02
%}
Mask_full = BET(iMag_mask,matrix_size_mask,voxel_size);
Mask_pad=padarray(Mask_full,[0,0,(matrix_size(3)-mask_red_num)],0,'post');

%This erodes the mask by a defined amount, this helps remove the skull
%effects on the background removal and thus the resulting QSM. Helps to
%remove streaking

images.MaskErode=2;
StEl=strel('cube',images.MaskErode); %try the cube structural element with n pixel erosion
Mask=imerode(Mask_pad,StEl);
images.Mask_red=Mask;

%Prepare mask based on magnitude thresholding
%Mask = genMask(iField,voxel_size);
 
mask_4D=repmat(Mask, 1,1,1, num_echo);

iFreq_raw=mask_4D.*PhaseData;
images.iFreq_raw=iFreq_raw;

%The mask on iField
Impose_Mask=abs(iField).*Mask;
images.Impose_Mask=Impose_Mask;

%The mask subtracted from iField
Neg_Mask=abs(iField).*(~Mask);
images.Neg_Mask=Neg_Mask;

%% Unwrapping
%{  
    Spatial unwrapping and temporal unwrapping
    unwrap each echo first and then combine the echoes together. This is done by linear fitting.
    Unwrapping is done to the phase and magnitude data
        Loop over different echoes and unwrap using the method of choice (this can be turned to a function)
        Combine using  CombinationEcho()

    There are several unwrapping techniques:
    1. SEGUE
    2. MEDI RG
    3. MEDI Laplacian
    4. Graph Cuts Subsample
    5. Graph Cuts
    6. Guillaume-Laplacian
    7. Quality-Guided-Vero

%}
% It appears that the issue of the negative multiplication is in the fit
% ppm - or so it seems. I will be using Guillames implementation to get
% around this. His implementation requires you to unwrap each echo first
% and then combine the echoes together. This is done by linear fitting.
% Since the wrapping is different for each echo it doesn't make sense to
% try combine them before unwrapping.. it makes me question the
% Fit_ppm_complex all together.


% For SEGUE
Inputs.Mask = Mask;
Inputs.Phase = PhaseData;
[PhaseData]=SEGUE(Inputs); images.UnwrapType='SEGUE';  


% For all other unwrapping techniques
% for jj=1:size(TE,2)    
%     display('Phase unwrapping')
% %     [UPhaseData]=unwrapPhase(MagnData(:,:,:,jj), iFreq_raw(:,:,:,jj), matrix_size); images.UnwrapType='MEDI RG';                               %MEDI RG
% %     [UPhaseData]=unwrapLaplacian(iFreq_raw(:,:,:,jj), matrix_size, voxel_size); images.UnwrapType='MEDI Laplacian';                            %MEDI Lap
% %     [UPhaseData]=unwrapping_gc(iFreq_raw(:,:,:,jj), iMag, voxel_size,2); images.UnwrapType='Graph Cuts Subsample';                               %Graph Cuts        
% %     [UPhaseData]=unwrapping_gc(iFreq_raw(:,:,:,jj), iMag, voxel_size); images.UnwrapType='Graph Cuts';                                         %Graph Cuts        
% %     [UPhaseData]=LaplacianUnwrap(PhaseData(:,:,:,jj),[voxel_size(1) voxel_size(2)],voxel_size(3)); images.UnwrapType='Guillaume-Laplacian';    %Gilbert
%     [UPhaseData]=qualityGuidedUnwrapping(PhaseData(:,:,:,jj),Mask);  images.UnwrapType='Quality-Guided-Vero';                                  %Veronique
% 
%     PhaseData(:,:,:,jj)=UPhaseData; 
%  end
% clear UPhaseData;

tic
disp('Temporal Unwrap Started')
PhaseData = AlexTemporalUnwrapping3(PhaseData, Mask);
images.TU='TUused';
toc
disp('Temporal Unwrap Finished')

% images.TU='TUnotused';

PhaseData=mask_4D.*PhaseData;
MagnData=mask_4D.*MagnData;
images.PhaseData=PhaseData;

disp('Combination of the echoes');
[iFreq, Resids]=CombinationEchoResiduals(MagnData,PhaseData,dataScan);
% [iFreq, Resids]=CombinationEcho(MagnData,PhaseData,dataScan);

iFreq=iFreq.*delta_TE;
images.iFreq=iFreq;
images.Residuals=Resids;



%% Background removal
%{
    Several methods of background removal exists. The list is presented
    below. More details in each section
        1. Projection onto Dipole Fields
        2. Laplacian Boundary Value
        3. sophisticated harmonic artifact reduction on phase data (SHARP)
        4. Regularized Enabled SHARP (RESHARP)
%}
N_std=1./iMag;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Background field removasl using Projection onto Dipole Fields

% RDF=iFreq;
% images.RDFtype='none';


RDF = PDF(iFreq, N_std, Mask ,matrix_size,voxel_size, B0_dir);
images.RDFtype = 'PDF';

%%%% NMR Biomed 2011;24(9):1129-36.
%%%% MRM 2010;63(1):194-206

% RDF = PDF_Gilbert(iFreq, N_std, Mask ,matrix_size,voxel_size, B0_dir);
% images.RDFtype = 'PDF-Gilbert';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Background field removal using Laplacian Boundary Value

% RDF = LBV(iFreq,Mask,matrix_size,voxel_size);
% images.RDFtype = 'LBV';

% Background field removal using Laplacian Boundary Value - Guillame
% implementation

% RDF = LBV_Gilbert(iFreq,Mask,matrix_size,voxel_size);
% images.RDFtype = 'LBVGilbert';


%%%% NMR Biomed 2014;27(3):312-319

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Background removal using sophisticated harmonic artifact reduction on
%phase data (SHARP)

% RDF = SHARP(iFreq, Mask, matrix_size, voxel_size);%, radius,threshold);
% images.RDFtype = 'SHARP';

%In SHARP there are a total of 6 parameters. The fifth is the radius 
%which is the radius of the spherical mean value operation and if it hasn't
%been specified, it defaults to: round(6/max(voxel_size)) * max(voxel_size)
%The 6th is a threshold parameter which is the threshold used in Truncated 
%SVD. If it hasn't been specified, it defaults to 0.00.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Background removal using Regularized Enabled SHARP (RESHARP)
% RDF = RESHARP(iFreq, Mask, matrix_size, voxel_size);% , radius,alpha);
% images.RDFtype = 'RESHARP';

%In RESHARP there are a total of 6 parameters. The fifth is the radius 
%which is the radius of the spherical mean value operation and if it hasn't
%been specified, it defaults to: round(6/max(voxel_size)) * max(voxel_size)
%The 6th is an alpha parameter which is the regularizaiton parameter used 
%in Tikhonov. If it hasn't been specified, it defaults to 0.01.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

images.RDF=RDF;

%% R2* and CSF Mask
% R2* map needed for ventricular CSF mask
R2s = arlo(TE, abs(iField));
images.R2s=R2s;
% Ventricular CSF mask for zero referencing 
%	Requirement:
%		R2s:	R2* map
Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
images.Mask_CSF=Mask_CSF;

%% Save stuff


save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
    voxel_size delta_TE CF B0_dir Mask_CSF;
 
%% QSM
%{
    QSM can be generated using several methods as described bellow
	� Morphology enabled dipole inversion with zero reference using CSF (MEDI+0)
	� iterative Susceptibility Weighted Imaging and Mapping (iSWIM)
	� iSWIM in veronique's implementations 
	� Truncated K-space Division (TKD) 
	� Truncated K-space Division (TKD) using zero referencing 
	� Truncated Singular Value Decomposition (TSVD)
	� Total Variation with Split Bregmann (TVSB)

%}

% Morphology enabled dipole inversion with zero reference using CSF (MEDI+0)
% images.QSMtype='MEDI';
% QSM = MEDI_L1('lambda',1000,'lambda_CSF',100);

% iterative Susceptibility Weighted Imaging and Mapping (iSWIM) 
% QSM = iSWIM(.15,0,'filename','RDF.mat');
% images.QSMtype='SWIM';

% iSWIM - Veroniques Implementation
% Define the forward and inverse filter in the Fourier domain with a threshold
% 
[ forwardFilter, inverseFilterReg ] = defineFilters( B0_dir,matrix_size,voxel_size );
% Initial QSM
susceptibilityMap=(ifftn(ifftshift(inverseFilterReg.*fftshift(fftn(RDF))))./(gamma*mainField*delta_TE).*10^6).*Mask;
% SWIM iterations
[QSMint,iterateInverse]=inverseProcess(real(susceptibilityMap), Mask, forwardFilter, matrix_size);

ValToMean=QSMint(logical(Mask_CSF));
QSM = QSMint - mean(ValToMean);
images.QSMtype='iSWIMvero';
% images.QSMtype='SWIMzeroREF';


% Truncated K-space Division (TKD) 
% QSM=TKD(.2,'filename','RDF.mat');
% images.QSMtype='TKD';
 
    
%QSM=Zero_ref_TKD(.2,'filename','RDF.mat');
%images.QSMtype='Zero-ref-TKD';


% Truncated Singular Value Decomposition (TSVD)
% QSM=TSVD(0.2,'filename','RDF.mat');
% images.QSMtype='TSVD';

% Total Variation with Split Bregmann (TVSB)
% QSM=TVSB('lambda',.0005,'filename','RDF.mat');
% images.QSMtype='TVSB';
% %


images.QSM=QSM;
images.PrettyQSM=QSM.*Mask;



%% Plot the data
Plotter2(images,FigSaveLocation);

%% Time the process
end_time=etime(clock, start_time);
images.end_time=end_time;

hours=end_time/3600;
minutes=60*(hours-floor(hours));
seconds=60*(minutes-floor(minutes));
if end_time < 60
    sprintf('The elapsed time is %.2f seconds', end_time)
elseif end_time < 3600
    sprintf('The elapsed time is %.0f minutes, and %.2f seconds', floor(minutes), seconds)
else 
    sprintf('The elapsed time is %.0f hours, %.0f minutes, and %.2f seconds', floor(hours), floor(minutes), seconds)
end
%% Save Data
save([DataSaveLocation, images.dataset '_Temp_Unwrapping_' num2str(images.TU) '_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype], 'images', '-v7.3');
