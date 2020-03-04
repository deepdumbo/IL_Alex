%**************************************************************************
% Implementation of the Laplacian phase unwrapping algorithm (Li et al.)
%**************************************************************************
function [Original_phase]=LaplacianUnwrap(Original_phase,AcqSpacing,SliceSpacing)

size_data=size(Original_phase);

%**************************************************************************
% Calculation of the k^2 matrix (we assume an isotropic inplane resolution)
%**************************************************************************
k2=zeros(size_data);
for k=1:size_data(1)
    for m=1:size_data(2)
        for n=1:size_data(3)
            k2(k,m,n)=(double(k)-double(floor(size_data(1)/2))-1.0)^2+(double(m)-double(floor(size_data(2)/2))-1.0)^2+((double(n)-double(floor(size_data(3)/2))-1.0).*(SliceSpacing./AcqSpacing(1)))^2;      
        end
    end
end

%**************************************************************************
% Equation 13 from the paper from Li et al.
%**************************************************************************
Original_phase=fftn(cos(Original_phase).*(ifftn((ifftshift(k2).*fftn(sin(Original_phase)))))-sin(Original_phase).*(ifftn((ifftshift(k2).*fftn(cos(Original_phase))))))./ifftshift(k2);

%**************************************************************************
% To prevent errors arising from the k=0 point
%**************************************************************************
Original_phase(isnan(Original_phase))=0; 
Original_phase(isinf(Original_phase))=0;

%**************************************************************************
% Back to the image domain
%**************************************************************************
%Original_phase=ifftn(Original_phase).*Mask;
Original_phase=ifftn(Original_phase);

end

