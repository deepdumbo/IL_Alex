%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWIM: streaking correction                                              %
% Véronique Fortier                                                       %
% 2015                                                                    %
% Based on: Tang, J., Liu, S., Neelavalli, J., Cheng, Y. C. N., Buch, S., %
% & Haacke, E. M. (2013). Improving susceptibility mapping using a        %
% threshold?based K?space/image domain iterative reconstruction approach. %
% Magnetic resonance in medicine, 69(5), 1396-1407.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ iterateQSM,iterateInverse ] = inverseProcess( susceptibilityMap, maskAirBone, forwardFilter,matrixSize)

    iterateInverse=0;
    thresholdInverse=0.004;
    thresholdFilter=0.1;

    % Define the mask for the undefined regions of the filter
    maskCone=forwardFilter;
    maskCone(find(abs(forwardFilter)<=thresholdFilter))=1;
    maskCone(find(abs(forwardFilter)>thresholdFilter))=0;
    
    % ___________________________________________________________________
    % Iterative inverse process
    while iterateInverse<=4     % I have to force the algorithm to stop after 4 iterations unless it diverges
%     while betaInverse>=thresholdInverse

        QSMRegionInterest=susceptibilityMap.*maskAirBone;
%         QSMRegionInterest=edge_preserving_avg(QSMRegionInterest,maskAirBone,4,voxelSize(1),voxelSize(2),voxelSize(3));
    
        FT_QSMRegionInterest=fftshift(fftn(QSMRegionInterest));

        Cone_FT_QSMRegionInterest=FT_QSMRegionInterest.*maskCone;

        FT_initialQSM=fftshift(fftn(susceptibilityMap));

        FT_initialQSM_noCone=FT_initialQSM.*(1-maskCone);
    
        FT_QSM_merged=FT_initialQSM_noCone+Cone_FT_QSMRegionInterest;
        
        if mod(matrixSize(3),2)==0
            FT_QSM_merged((matrixSize(1)/2+1),(matrixSize(2)/2+1),(matrixSize(3)/2+1))=0;
        end
        iterateQSM=ifftn(ifftshift(FT_QSM_merged));
    
        somme=sum(sum(sum((real(susceptibilityMap)-real(iterateQSM)).^2)));
        betaInverse=sqrt(somme/(matrixSize(1).*matrixSize(2).*matrixSize(3)));
        

        susceptibilityMap=iterateQSM;
               
        iterateInverse=iterateInverse+1;
   
        if betaInverse<=thresholdInverse
            break
        end
    end


end

