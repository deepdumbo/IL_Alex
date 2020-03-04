%% Load file

% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end
%% Initialization 

gamma=267.513E6;
CF=images.CF;
voxel_size=images.voxel_size;
mainField=images.mainField;
TE=images.TE;
delta_TE=images.delta_TE;
B0_dir=images.B0_dir;
matrix_size=images.matrix_size;
Mask=images.Mask_red;
iFreq=images.iFreq;
iMag=images.iMag;
N_std=1./iMag;
RDF=images.RDF;
iFreq_raw=RDF; %I never save iFreq_raw.. so this could be an issue?
R2s=images.R2s;
Mask_CSF=images.Mask_CSF;

%% Background removal

% RDF = PDF(iFreq, N_std, Mask ,matrix_size,voxel_size, B0_dir);
% images.RDFtype = 'PDF';

% RDF = PDF_Gilbert(iFreq, N_std, Mask ,matrix_size,voxel_size, B0_dir);
% images.RDFtype = 'PDF-Gilbert';

% RDF = LBV_Gilbert(iFreq,Mask,matrix_size,voxel_size);
% images.RDFtype = 'LBVGilbert';

% RDF = LBV(iFreq,Mask,matrix_size,voxel_size);
% images.RDFtype = 'LBV';



images.RDF=RDF;
%% Save
save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
    voxel_size delta_TE CF B0_dir Mask_CSF;
%% Dipole inversion

% images.QSMtype='MEDI';
% QSM = MEDI_L1('lambda',1000,'lambda_CSF',100);

% QSM = iSWIM(.15,0,'filename','RDF.mat');
% images.QSMtype='SWIM';

% [ forwardFilter, inverseFilterReg ] = defineFilters( B0_dir,matrix_size,voxel_size );
% Initial QSM
susceptibilityMap=(ifftn(ifftshift(inverseFilterReg.*fftshift(fftn(RDF))))./(gamma*mainField*delta_TE).*10^6).*Mask;
% SWIM iterations
[QSM,iterateInverse]=inverseProcess(real(susceptibilityMap), Mask, forwardFilter, matrix_size);
% images.QSMtype='iSWIMvero';
QSMint=QSM;
ValToMean=QSMint(logical(Mask_CSF));
QSM = QSMint - mean(ValToMean);
images.QSMtype='SWIMzeroREF';

% DidItWork=mean(QSM(logical(Mask_CSF)));

% Truncated K-space Division (TKD) 
% QSM=TKD(.2,'filename','RDF.mat');
% images.QSMtype='TKD';

% Truncated K-space Division (TKD) using zero referencing
% QSMint=TKD(.2,'filename','RDF.mat');
% ValToMean=QSMint(logical(Mask_CSF));
% QSM = QSMint - mean(ValToMean);
% images.QSMtype='TKDzeroREF';

% QSM=TSVD(0.2,'filename','RDF.mat');
% images.QSMtype='TSVD';

% QSM=TVSB('lambda',.0005,'filename','RDF.mat');
% images.QSMtype='TVSB';

images.QSM=QSM;
images.PrettyQSM=QSM.*Mask;

%% Plot the data
Plotter2(images);
%% Save Data
save(['./Saved_MAT_files/' images.dataset '_Mask_erosion_' num2str(images.MaskErode) '_' images.UnwrapType '_' images.RDFtype '_' images.QSMtype], 'images', '-v7.3');