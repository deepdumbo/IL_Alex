%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal unwrapping                                                      %
% Véronique Fortier                                                       %
% 2017                                                                    %
% Input: 
%   -phaseModelUW: 4D phase matrix, 4th dimension corresponding to echo
%   time. All echoes must be spatially unwrapped first.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ phaseModelUW_t ] = temporalUnwrapping( phaseModelUW )


phaseModelUW_t=phaseModelUW;


for i=1:size(phaseModelUW,4)-1
    diff=phaseModelUW_t(:,:,:,i+1)-phaseModelUW_t(:,:,:,i);
    result=phaseModelUW_t(:,:,:,i+1);
    
    while isempty(find(diff>=(2*pi)))==0 
        result(find(diff>=(2*pi)))=result(find(diff>=(2*pi)))-2*pi;
        diff=result-phaseModelUW_t(:,:,:,i);
    end
    
    while isempty(find(diff<=-(2*pi)))==0 
        result(find(diff<=-(2*pi)))=result(find(diff<=-(2*pi)))+2*pi;
        diff=result-phaseModelUW_t(:,:,:,i);
    end

    phaseModelUW_t(:,:,:,i+1)=result;
end


end

