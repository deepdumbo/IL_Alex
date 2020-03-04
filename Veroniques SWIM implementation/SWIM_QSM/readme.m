% SWIM QSM


%% STEP 1 : green function

% Define the forward and inverse filter in the Fourier domain with a threshold
[ forwardFilter, inverseFilterReg ] = defineFilters( B0_dir,matrix_size,voxel_size );


% Initial QSM
susceptibilityMap=(ifftn(ifftshift(inverseFilterReg.*fftshift(fftn(RDF))))./(gamma*mainField*delta_TE).*10^6).*Mask;


% SWIM iterations

[iterateQSM,iterateInverse]=inverseProcess(real(susceptibilityMap), Mask, forwardFilter, matrix_size);