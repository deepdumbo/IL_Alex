%**************************************************************************
% Slice by slice brain masking function
%**************************************************************************
function [Mask]=AutoMaskRec(Complex_data)

size_mask=size(Complex_data);
Mask_magn=zeros(size_mask);
for k=1:size_mask(3)
    maxim=max(max(abs(Complex_data(:,:,k))));
    temp=abs(Complex_data(:,:,k))./maxim;
    Magn_thresh=graythresh(temp);
    Mask_magn(:,:,k)=im2bw(temp,Magn_thresh);    
end

Mask=Mask_magn;

%**************************************************************************
% Identification of the main label (which should be the brain)
%**************************************************************************

%**************************************************************************
%Light erosion
%**************************************************************************
SE1=ones(2,2,2);
Mask=imerode(Mask,SE1);

for k=1:size(Mask,3)
    %**********************************************************************
    %Identification of all connected regions
    %**********************************************************************
    [Label_mask num]=bwlabel(Mask(:,:,k));
    
    %**********************************************************************
    %Identification of the largest region
    %**********************************************************************
    Size_label=zeros(num,1);
    
    for m=1:num
        Label_mask2=Label_mask(:);
        Label_mask2(Label_mask2~=m)=[];
        Size_label(m)=length(Label_mask2);
    end
    
    if isempty(Size_label)==0
        [num Label]=max(Size_label);
        Label_mask(Label_mask~=Label)=0;
        Label_mask(Label_mask==Label)=1;
    else
        Label_mask=zeros(size_mask(1),size_mask(2));
    end
    %Mask(:,:,k)=imfill(Mask(:,:,k));
    Mask(:,:,k)=imfill(Label_mask);
end

%**************************************************************************
%Morphological closing operation
%**************************************************************************
r=3;
[x,y,z] = meshgrid(-r:r,-r:r,-r:r);
SE2 = (x/r).^2 + (y/r).^2 + (z/r).^2 <= 1;
Mask=imclose(Mask,SE2);

%**************************************************************************
% Light dilatation
%**************************************************************************
%Mask=imdilate(Mask,SE1);


end



