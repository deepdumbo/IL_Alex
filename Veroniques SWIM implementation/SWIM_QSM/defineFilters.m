%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dipole response                                                         %
% Véronique Fortier                                                       %
% 2015                                                                    %
% Based on: Haacke, E. M., Liu, S., Buch, S., Zheng, W., Wu, D., & Ye, Y. %
% (2015). Quantitative susceptibility mapping: current status and future  %
% directions. Magnetic resonance imaging, 33(1), 1-25.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ forwardFilter, inverseFilterReg ] = defineFilters( B0_dir,matrixSize,voxelSize )

clear kZ kX ky kx kz

thresholdFilter=0.1;


fov=matrixSize(1)*voxelSize(1);     %  in x and y plane
fovZ=matrixSize(3)*voxelSize(3);

kX=-(1/voxelSize(1))/2:(1/fov):(1/voxelSize(1))/2-1/fov;
kZ(1,1,:)=-(1/voxelSize(3))/2:1/fovZ:(1/voxelSize(3))/2-1/fovZ;

% kx and ky are constructed from identical pixel size in each direction, so
%   in effect, voxelSize(1) and voxelSize(2) are inherently assumed to be
%   equal here
ky=repmat(kX',1,matrixSize(1),matrixSize(3));
kx=repmat(kX,matrixSize(1),1,matrixSize(3));
kz=repmat(kZ,[matrixSize(1),matrixSize(1),1]);

inverseFilterReg=zeros(matrixSize(1),matrixSize(1),matrixSize(3));
% inverseFilterReg2=zeros(matrixSize(1),matrixSize(1),matrixSize(3));
forwardFilter=zeros(matrixSize(1),matrixSize(1),matrixSize(3));

denominator=kx.^2+ky.^2+kz.^2;

% Note that the definition of the foward filter "defines" the direction of
% the applied magnetic field (k)
if B0_dir==[0 1 0]
    k=ky;
    other1=kx;
    other2=kz;
elseif B0_dir==[0 -1 0]
    k=ky;
    other1=kx;
    other2=kz;
elseif B0_dir==[1 0 0]
    k=kx;
    other1=kz;
    other2=ky;
elseif B0_dir==[-1 0 0]
    k=kx;
    other1=kz;
    other2=ky;
elseif B0_dir==[0 0 -1]
    k=kz;
    other1=kx;
    other2=ky;
else
    k=kz;
    other1=kx;
    other2=ky;
end

forwardFilter(find(denominator==0))=0;
indexForward=find(denominator~=0);
forwardFilter(indexForward)=(1/3)-((k(indexForward).^2)./denominator(indexForward));

index=find(abs((1/3)-((k.^2)./denominator))>thresholdFilter);
inverseFilterReg(index)=1./((1/3)-((k(index).^2)./denominator(index)));
% inverseFilterReg2(index)=1./((1/3)-((k(index).^2)./denominator(index)));

% Identify regions to regularize
indexReg=find(abs((1/3)-((k.^2)./denominator))<=thresholdFilter);



%% Method in the Haacke review 2015 for regularization, bring slowly to 0
% inverseFilterReg(indexReg)=thresholdFilter^(-3).*sign((1/3)-((k(indexReg).^2)./denominator(indexReg))).*((1/3)-((k(indexReg).^2)./denominator(indexReg))).^2;
% inverseFilterReg(indexReg)=0;


%% Method Haacke et al, 2010
% First part: set g^-1 to 1/threshold, keep the sign
indexNegative=find((((1/3)-((k.^2)./denominator))>-thresholdFilter)&(((1/3)-((k.^2)./denominator))<0));
indexPositive=find((((1/3)-((k.^2)./denominator))<thresholdFilter)&(((1/3)-((k.^2)./denominator))>=0));

inverseFilterReg(indexNegative)=-1/thresholdFilter;
inverseFilterReg(indexPositive)=1/thresholdFilter;
inverseFilterReg(isnan(inverseFilterReg)) =0;
inverseFilterReg(isinf(inverseFilterReg)) =0;

% Second part: Bring slowly to zero by multiplying g^-1 by alpha^2
Z0 = sqrt( (other1.^2+other2.^2)/2);
alpha = abs(abs(k)-Z0);
positive_Za = sqrt((1-3*thresholdFilter)*(other1.^2+other2.^2)/(2+3*thresholdFilter));
negative_Za = sqrt((1+3*thresholdFilter)*(other1.^2+other2.^2)/(2-3*thresholdFilter));

negative_region = (forwardFilter<0).*alpha./abs(negative_Za - Z0);
positive_region = (forwardFilter>0).*alpha./abs(positive_Za - Z0);
alpha = negative_region+positive_region;
alpha(alpha>1) = 1;

inverseFilterReg=inverseFilterReg.*alpha.^2;
inverseFilterReg(isnan(inverseFilterReg)) = 0;    
inverseFilterReg(isinf(inverseFilterReg)) = 0;

indexCenter=find(denominator==0);
if isempty(indexCenter)~=1
    inverseFilterReg(indexCenter)=0;
end


%% from guillaume code

% forwardFilter = dipole_kernelKS(matrixSize, voxelSize, 1,'kspace');
% 
% thresholdFilter = 1/thresholdFilter;
% inverseFilterReg = 1/forwardFilter;
% inverseFilterReg((inverseFilterReg>thresholdFilter)) = thresholdFilter;
% inverseFilterReg((1/forwardFilter<-thresholdFilter)) = -thresholdFilter;
% 
% 


end

