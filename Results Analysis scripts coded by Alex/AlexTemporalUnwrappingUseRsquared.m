
%{
    last update 2020/03/06, Hossein
    *Does not work for now*

    a Function that does temporal phase unrwapping.
    input:  images with spatially unwrapped phase
            Mask
            
    Output: images unwrapped in time
%}

% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end


% This code works as it is suppose to, however the logic behind it is
% flawed slightly resulting in it not being very robust.

%The issue is that by modifying the fit, the r^2 doesn't get significantly
%better when it should.
function [ TempUnwrappedPhase] = AlexTemporalUnwrappingUseRsquared(SpatUnwrappedPhase,TE, Mask)

iter=1;
% Step 0:
% Get the median value of the entire 3D volume for each echo.
numechoes = size(SpatUnwrappedPhase,4);
numslices = size(SpatUnwrappedPhase,3);
TheMedians=zeros(1,numechoes);

h=45;
% for h=1:numslices
    for i=1:numechoes
        TheMedians(i) = median(nonzeros(SpatUnwrappedPhase(:,:,h,i)));
    end
    disp(TheMedians)
% Step 1:
% Plot these mean values against their echo time and get a linear fit. 
    while iter<10
        x = TE;
        y = TheMedians;
        [TheFit, S] = polyfit(x,y,1);
        % Record the R2, or some other value to measure the goodness of fit.
        f = polyval(TheFit,x);
        plot(x,y,'o',x,f,'-')
        Rsquared = 1 - (S.normr/norm(y - mean(y)))^2;
        disp(Rsquared)


        [~,MaxLoc] = max(TheMedians);
        PhaseDataMaxRed = SpatUnwrappedPhase;
        PhaseDataMaxRed(:,:,h,MaxLoc) = SpatUnwrappedPhase(:,:,h,MaxLoc)-(2*pi);
        PhaseDataMaxRed = PhaseDataMaxRed.*Mask;
        TheMediansMax = zeros(1,numechoes);
        for i=1:numechoes
            TheMediansMax(i) = median(nonzeros(PhaseDataMaxRed(:,:,h,i)));
        end
        disp(TheMediansMax)
        y_max = TheMediansMax;
        [TheFit_max, S_max] = polyfit(x,y_max,1);
        % Record the R2, or some other value to measure the goodness of fit.
        f = polyval(TheFit_max,x);
        plot(x,y_max,'o',x,f,'-')
        Rsquared_max = 1 - (S_max.normr/norm(y_max - mean(y_max)))^2;
        disp(Rsquared_max)
        
        
        [~,MinLoc] = min(TheMedians);
        PhaseDataMinInc = SpatUnwrappedPhase;
        PhaseDataMinInc(:,:,h,MinLoc) = SpatUnwrappedPhase(:,:,h,MinLoc)+(2*pi);
        PhaseDataMinInc = PhaseDataMinInc.*Mask;
        TheMediansMin = zeros(1,numechoes);
        for i=1:numechoes
            TheMediansMin(i) = median(nonzeros(PhaseDataMinInc(:,:,h,i)));
        end
        y_min = TheMediansMin;
        disp(TheMediansMin)
        [TheFit_min, S_min] = polyfit(x,y_min,1);
        % Record the R2, or some other value to measure the goodness of fit.
        f = polyval(TheFit_min,x);
        plot(x,y_min,'o',x,f,'-')
        Rsquared_min = 1 - (S_min.normr/norm(y_min - mean(y_min)))^2;
        disp(Rsquared_min)
        
        if Rsquared_max > Rsquared && Rsquared_max > Rsquared_min
            SpatUnwrappedPhase=PhaseDataMaxRed;
        elseif Rsquared_min > Rsquared && Rsquared_min > Rsquared_max
            SpatUnwrappedPhase=PhaseDataMinInc;
        else
            breaks
        end
        iter=iter+1;
    end
    

% end
disp(TheMedians)
% TempUnwrappedPhase=SpatUnwrappedPhase;



end
