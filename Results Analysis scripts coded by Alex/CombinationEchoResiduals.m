%**************************************************************************
% Implementation of the combination of phase images from different echo using a
% weighted linear regression (Gilbert et al.)
%**************************************************************************
function [MapdB, Residuals]=CombinationEchoResiduals(MagnData,PhaseData,info)

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

%**************************************************************************
% We smooth the phase offset using a 4th order polynomial filter
%**************************************************************************
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
NoiseMap=sum(1./sigma.^2,4)./((sum(1./sigma.^2,4).*sum(TEm.^2./sigma.^2,4))-sum(TEm./sigma.^2,4).^2);
NoiseMap=sqrt(NoiseMap);

% Apparently I am not suppose to find the closest distance, but instead the
% vertical distance. This is the closest distance.
%
% %**************************************************************************
% % Alex's edits to get a map of residuals - the residuals that describe the
% % goodness of fit
% %**************************************************************************
% disp('Starting residual calculation...')
% tic
% D=zeros(1,5);
% Residuals=zeros(size(PhaseData,1),size(PhaseData,2),size(PhaseData,3));
% % Guaranteed, there is a faster way to do this - matrix wise
% for g=1:size(PhaseData,1)
%     for h=1:size(PhaseData,2)
%         for i=1:size(PhaseData,3)
%             for k=1:length(info.TE)
%                 %What is the point you are trying to compare?
%                 x1=info.TE(k);
%                 y1=PhaseData(g,h,i,k);                
%                 OrigSlope = MapdB(g,h,i); %The slope of the perpendicular line is the negative inverse of the original slope               
%                 
%                 NewLineInt = y1+(x1/OrigSlope);
%                 x2 = (NewLineInt)/(OrigSlope+(1/OrigSlope));
%                 y2 = x2*OrigSlope;
%                 
%                 %Find the distance using pythagorean theorem
%                 D(k)=sqrt((x2-x1)^2+(y2-y1)^2);
%                                                          
%             end
%             %Sum the distances for all echoes in one voxel
%             Residuals(g,h,i) = sum(D);
%         end
%     end
% end
% toc


%**************************************************************************
% Alex's edits to get a map of residuals - the residuals that describe the
% goodness of fit
%**************************************************************************

disp('Starting residual calculation...')
tic
D=zeros(1,5);
Residuals=zeros(size(PhaseData,1),size(PhaseData,2),size(PhaseData,3));
% Guaranteed, there is a faster way to do this - matrix wise
for g=1:size(PhaseData,1)
    for h=1:size(PhaseData,2)
        for i=1:size(PhaseData,3)
            for k=1:length(info.TE)
                x1=info.TE(k);
                y1=PhaseData(g,h,i,k);   %Phi0 has already been subtracted - see line 33             
                y2 = x1*MapdB(g,h,i); %Find point on the model with same TE as the theoretical value                                              
                %Find the residual
                D(k)=(y2-y1)^2;                                                        
            end
            %Sum the distances for all echoes in one voxel
            Residuals(g,h,i) = sqrt(sum(D));
        end
    end
end
toc





%End for the function
end


