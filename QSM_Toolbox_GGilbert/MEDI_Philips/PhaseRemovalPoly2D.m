%**************************************************************************
% Polynomial filtering (3rd order)
%**************************************************************************

function OutputImage=PhaseRemovalPoly2D(OriginalImage)

Dim=size(OriginalImage);

%**************************************************************************
%Zero-th order term
%**************************************************************************
Const=ones(Dim);

%**************************************************************************
%First order
%**************************************************************************
[X Y]=meshgrid(-floor(Dim(2)/2):ceil(Dim(2)/2)-1,-floor(Dim(1)/2):ceil(Dim(1)/2)-1);

%**************************************************************************
%2nd order
%**************************************************************************
X2=X.*X;
Y2=Y.*Y;
XY=X.*Y;


%**************************************************************************
%3rd order
%**************************************************************************
X3=X.*X.*X;
Y3=Y.*Y.*Y;
X2Y=X.*X.*Y;
Y2X=Y.*Y.*X;

%**************************************************************************
% 4th order
%**************************************************************************
X4=X.*X.*X.*X;
Y4=Y.*Y.*Y.*Y;
X3Y=X.*X.*X.*Y;
X2Y2=X.*X.*Y.*Y;
Y3X=X.*Y.*Y.*Y;

%**************************************************************************
%Mask
%**************************************************************************
Mask=zeros(Dim);
Mask(OriginalImage~=0)=1;
Mask=Mask(:);

B=OriginalImage(:).*Mask;
B(Mask==0)=[];

%First order
temp=Const(:);
temp(Mask==0)=[];
A(:,1)=temp;
temp=X(:);
temp(Mask==0)=[];
A(:,2)=temp;
temp=Y(:);
temp(Mask==0)=[];
A(:,3)=temp;

%Second order
temp=X2(:);
temp(Mask==0)=[];
A(:,4)=temp;
temp=Y2(:);
temp(Mask==0)=[];
A(:,5)=temp;
temp=XY(:);
temp(Mask==0)=[];
A(:,6)=temp;

%Third order
temp=X3(:);
temp(Mask==0)=[];
A(:,7)=temp;
temp=Y3(:);
temp(Mask==0)=[];
A(:,8)=temp;
temp=X2Y(:);
temp(Mask==0)=[];
A(:,9)=temp;
temp=Y2X(:);
temp(Mask==0)=[];
A(:,10)=temp;

%Fourth order
temp=X4(:);
temp(Mask==0)=[];
A(:,11)=temp;
temp=Y4(:);
temp(Mask==0)=[];
A(:,12)=temp;
temp=X3Y(:);
temp(Mask==0)=[];
A(:,13)=temp;
temp=X2Y2(:);
temp(Mask==0)=[];
A(:,14)=temp;
temp=Y3X(:);
temp(Mask==0)=[];
A(:,15)=temp;

%**************************************************************************
%We find the solution to this linear algebra problem
%**************************************************************************
Coeff=A\B;

FilteredImage=(Coeff(1)*Const+Coeff(2)*X+Coeff(3)*Y+Coeff(4)*X2+Coeff(5)*Y2+Coeff(6)*XY+...
    Coeff(7)*X3+Coeff(8)*Y3+Coeff(9)*X2Y+Coeff(10)*Y2X+Coeff(11)*X4+Coeff(12)*Y4+Coeff(13)*X3Y+Coeff(14)*X2Y2+Coeff(15)*Y3X);

OutputImage=FilteredImage;

end