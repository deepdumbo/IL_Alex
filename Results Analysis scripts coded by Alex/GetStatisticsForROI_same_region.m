%{
last update 2020/03/04, Hossein

This script allows you to draw a ROI on your image and immediately learn some info on the ROI. 
Plot either the pretty QSM or ifield 
Ask the user if they want to use a previous ROI or draw a new one
	Depending on the answer, apply the mask and report the following:
	Numberofpixels, minGL, maxGL stdGL, meanGL
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just need to be sure that the correct image is loaded

load('C:\Users\alexe\OneDrive - McGill University\Alex MSc Thesis Research\QSM\Acquired data Nov 8 2019\Saved_MAT_files\Nov 18 - Investigating halfway start\Axial_Mask_erosion_4_Graph Cut_PDF_TKD.mat')


% Or choose the file yourself:

% [file,path] = uigetfile('*.mat');
% if isequal(file,0)
%    disp('User selected Cancel');
%    return
% else
%    load(fullfile(path,file));
% end

images.slice=65;


% Save the image being examined
% <HJ, Question: can this be turned into an if statement>
QSMimage = images.PrettyQSM(:,:,images.slice);
% QSMimage=abs(images.iField(:,:,images.slice));
 

% Plot the image
imagesc(QSMimage,[-0.3,0.3]);
axis image
title(sprintf('%s RDF, %s QSM, slice %d', images.RDFtype, images.QSMtype, images.slice))
set(gcf, 'Position', get(0,'Screensize')) %Maximize the image 
colorbar
impixelinfo
colormap(bone)


status=1;
counter=0;


answer = questdlg('Use previous ROI or draw a new one?', ... %The message in the dialog box
        'Previous or new?', ...                    %This is the title of the dialog box
        'Use previous', 'Draw new', 'Done');    %The options for the dialog box, not sure why Done needs to be there twice but it does..
    
switch answer %To let the while loop repeat
    case 'Use previous'
        
        counter=1;
        hFH = drawfreehand('Position', ROI_loc);
        binaryImage = hFH.createMask(); %creates a binary mask for your image
        %Get some information about the ROI
        numberOfPixels1 = sum(binaryImage(:)); %Count the number of pixels in your image


        blackMaskedImage = QSMimage; %Save it as something else so it doesn't get overwritten
        blackMaskedImage(~binaryImage) = 0; %Set all values outside of ROI to zero


        temp=cast(blackMaskedImage(binaryImage),'single'); %Save your ROI

        % Calculate some basic stats now
        meanGL = mean(temp);
        stdGL = std(temp);
        maxGL = max(temp);
        minGL = min(temp);

        %Display your results in a text box now
        message = sprintf('ROI # = %d\nMean value within drawn area = %.3f\nStandard Deviation = %.3f\nMax Value = %.2f\nMin Value = %.2f\nNumber of pixels = %d', ... 
            counter,meanGL,stdGL, maxGL, minGL, numberOfPixels1);
        set(msgbox(message),'Position',[440 462.5 150 100]); %[^right ^up ^wider ^taller] 
        return;
        
    case 'Draw new'
        while status == 1
        
            answer = questdlg('Press Done when finished, or Draw to draw an ROI', ... %The message in the dialog box
                'Continue?', ...                    %This is the title of the dialog box
                'Draw', 'Done', 'Done');    %The options for the dialog box, not sure why Done needs to be there twice but it does..

            switch answer %To let the while loop repeat
                case 'Draw'
                    status=1;
                    counter=counter+1;
                case 'Done'
                    status=0;
                    return; %Jump out of the while loop and terminate the program
            end


            % Draw your ROI
            %hFH = drawfreehand('Multiclick', true);
            hFH = drawfreehand();
            binaryImage = hFH.createMask(); %creates a binary mask for your image
            ROI_loc=hFH.Position;
            %Get some information about the ROI
            numberOfPixels1 = sum(binaryImage(:)); %Count the number of pixels in your image


            blackMaskedImage = QSMimage; %Save it as something else so it doesn't get overwritten
            blackMaskedImage(~binaryImage) = 0; %Set all values outside of ROI to zero


            temp=cast(blackMaskedImage(binaryImage),'single'); %Save your ROI

            % Calculate some basic stats now
            meanGL = mean(temp);
            stdGL = std(temp);
            maxGL = max(temp);
            minGL = min(temp);

            %Display your results in a text box now
            message = sprintf('ROI # = %d\nMean value within drawn area = %.3f\nStandard Deviation = %.3f\nMax Value = %.2f\nMin Value = %.2f\nNumber of pixels = %d', ... 
                counter,meanGL,stdGL, maxGL, minGL, numberOfPixels1);
            set(msgbox(message),'Position',[440 462.5 150 100]); %[^right ^up ^wider ^taller] 
        end

end

