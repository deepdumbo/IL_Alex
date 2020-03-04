%**************************************************************************
% Implementation of the combination of phase images from different echo using a
% weighted linear regression (Gilbert et al.)
%**************************************************************************
function [MapdB Residual]=CombinationEcho(MagnData,PhaseData,info)

 %**************************************************************
 % Initialization of some variables
 %**************************************************************
size_data=size(MagnData);

%**************************************************************************
% We estimate the phase offset from the first two echoes (much faster than
% a full regression
%**************************************************************************
Phi0=squeeze(PhaseData(:,:,:,1))-(squeeze(PhaseData(:,:,:,2)-PhaseData(:,:,:,1))/(info.TE(2)-info.TE(1)))*info.TE(1);
Phi0(isnan(Phi0))=0;
Phi0(isinf(Phi0))=0;

% %**************************************************************************
% % We smooth the phase offset using a 4th order polynomial filter
% %**************************************************************************
for k=1:size_data(3)
    Phi0(:,:,k)=PhaseRemovalPoly2D(Phi0(:,:,k));
end

%**************************************************************************
% Weighted linear regression 
%**************************************************************************
TEm=zeros(size_data(1),size_data(2),size_data(3),length(info.TE));
for k=1:length(info.TE)
    TEm(:,:,:,k)=info.TE(k);
    PhaseData(:,:,:,k)=PhaseData(:,:,:,k)-Phi0;
end
MapdB=sum(MagnData.^2.*PhaseData.*TEm,4)./sum(MagnData.^2.*TEm.^2,4);

%From Kressler et al.
%MapdB=((sum(MagnData.^2,4).*sum(MagnData.^2.*PhaseData.*TEm,4))-(sum(MagnData.^2.*TEm,4).*sum(MagnData.^2.*PhaseData,4)))./((sum(MagnData.^2,4).*sum(MagnData.^2.*TEm.^2,4))-(sum(MagnData.^2.*TEm,4).^2));

MapdB(isnan(MapdB))=0;
MapdB(isinf(MapdB))=0;

%**************************************************************************
% We identify the voxels that shown a unreliable phase by calculation of the
% field map standard deviation
%**************************************************************************

%From Kressler et al., assuming sigma=1;
sigma=1./MagnData;
Residual=sum(1./sigma.^2,4)./((sum(1./sigma.^2,4).*sum(TEm.^2./sigma.^2,4))-sum(TEm./sigma.^2,4).^2);
Residual=sqrt(Residual);


end


