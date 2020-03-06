%{
    last update 2020/03/06, Hossein
    *Does not work for now*

    a Function that does temporal phase unrwapping.
    input:  images with spatially unwrapped phase
            TE
            Mask
    Output: images unwrapped in time
%}
[file,path] = uigetfile('*.mat');
if isequal(file,0)
   disp('User selected Cancel');
   return
else
   load(fullfile(path,file));
end

function [ TempUnwrappedPhase] = AlexTemporalUnwrapping(SpatUnwrappedPhase,TE, Mask)


% Step 0:
% Get the mean value of the entire 3D volume for each echo.
TheMeans=zeros(1,size(SpatUnwrappedPhase,4));
numechoes = size(SpatUnwrappedPhase,4);

for i=1:numechoes
    TheMeans(i)=mean(mean(mean(SpatUnwrappedPhase(:,:,:,i))));
end

% Step 1:
% Plot these mean values against their echo time and get a linear fit. 

x = TE;
y = TheMeans;
[TheFit, S] = polyfit(x,y,1);
% Record the R2, or some other value to measure the goodness of fit.
f = polyval(TheFit,x);
plot(x,y,'o',x,f,'-')
RsquaredOne = 1 - (S.normr/norm(y - mean(y)))^2;

[MaxVal1,MaxLoc1] = max(TheMeans);
PhaseDataRed = SpatUnwrappedPhase;
PhaseDataRed(:,:,:,MaxLoc1) = SpatUnwrappedPhase(:,:,:,MaxLoc1)-(2*pi);
PhaseDataRed=PhaseDataRed.*Mask;
TheMeansMaxRed=zeros(1,size(SpatUnwrappedPhase,4));
for i=1:numechoes
    TheMeansMaxRed(i)=mean(mean(mean(PhaseDataRed(:,:,:,i))));
end

y = TheMeansMaxRed;
[TheFitOne, S] = polyfit(x,y,1);
f = polyval(TheFitOne,x);
plot(x,y,'o',x,f,'-')
RsquaredOneAdjusted = 1 - (S.normr/norm(y - mean(y)))^2;

% x = TE;
% y = TheMeans;
% [TheFit, S] = polyfit(x,y,1);
% % Record the R2, or some other value to measure the goodness of fit.
% f = polyval(TheFit,x);
% plot(x,y,'o',x,f,'-')
% RsquaredMaxRed = 1 - (S.normr/norm(y - mean(y)))^2;
% 
% x = TE;
% y = TheMeans;
% [TheFit, S] = polyfit(x,y,1);
% % Record the R2, or some other value to measure the goodness of fit.
% f = polyval(TheFit,x);
% plot(x,y,'o',x,f,'-')
% RsquaredMinInc = 1 - (S.normr/norm(y - mean(y)))^2;
end
