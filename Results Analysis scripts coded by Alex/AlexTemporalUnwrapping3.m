%{
    last update 2020/03/06, Hossein
    *Does not work for now*

    a Function that does temporal phase unrwapping.
    input:  images with spatially unwrapped phase
            Mask
            
    Output: images unwrapped in time
%}

function [ TempUnwrappedPhase] = AlexTemporalUnwrapping3(SpatUnwrappedPhase, Mask)


numechoes = size(SpatUnwrappedPhase,4);
TheMedians=zeros(1,numechoes);
numslices = size(SpatUnwrappedPhase,3);
diff = zeros(1,numechoes-1);
% h=45;

for h=1:numslices
    for i=1:numechoes
        TheMedians(i) = median(nonzeros(SpatUnwrappedPhase(:,:,h,i)));
    end


    for j = 1:(numechoes-1)
        diff(j)=TheMedians(j+1)-TheMedians(j);
    end
 

    for jj=1:(numechoes-1)
        if diff(jj)>(1.2*pi)
            SpatUnwrappedPhase(:,:,h,jj+1)=(SpatUnwrappedPhase(:,:,h,jj+1)-(2*pi)).*Mask(:,:,h);
            for i=1:numechoes
                TheMedians(i) = median(nonzeros(SpatUnwrappedPhase(:,:,h,i)));
            end
            for j = 1:(numechoes-1)
                diff(j)=TheMedians(j+1)-TheMedians(j);
            end
        elseif diff(jj)<(-1.2*pi)
            SpatUnwrappedPhase(:,:,h,jj+1)=(SpatUnwrappedPhase(:,:,h,jj+1)+(2*pi)).*Mask(:,:,h);
            for i=1:numechoes
                TheMedians(i) = median(nonzeros(SpatUnwrappedPhase(:,:,h,i)));
            end
            for j = 1:(numechoes-1)
                diff(j)=TheMedians(j+1)-TheMedians(j);
            end
        end
    end
end

TempUnwrappedPhase=SpatUnwrappedPhase;

end