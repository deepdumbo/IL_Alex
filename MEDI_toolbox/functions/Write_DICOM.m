function Write_DICOM(M,files,outdir,opts)
%WRITE_DICOM Summary of this function goes here
%   Detailed explanation goes here


defopts.SeriesDescription = 'QSM';
defopts.SeriesInstanceUID = [];
defopts.SeriesNumber = [];
defopts.Convert = @(x) convert(x);
defopts.Window = 500;
defopts.Level = 0;
% defopts.EchoNumber = [];
% defopts.EchoTime = 0.0;
        
if nargin<4; opts=struct; end
deffields=fieldnames(defopts);
for i=1:length(deffields)
    if ~isfield(opts, deffields{i})
        opts.(deffields{i})=defopts.(deffields{i});
    end
end

if isfield(files,'M')
    filenames=files.M;
elseif isfield(files,'R')
    filenames=files.R;
else
    error('No filenames (M nor R) found.');
end
    
flag_signed=min(M(:))<0;

if size(M,3) ~= size(filenames,1)
    error([num2str(size(filenames,1)) ' filenames given for ' num2str(size(M,3)) ' slices.']);
end

if isempty(opts.SeriesInstanceUID)
   opts.SeriesInstanceUID=dicomuid;
end
progress='';

mkdir(outdir);

warning('off','images:dicom_add_attr:invalidAttribChar');

for slice=1:size(M,3)
    info = dicominfo(filenames{slice,end});
    im = M(:,:,slice);
    if (isfield(info, 'SmallestImagePixelValue'))
        info.SmallestImagePixelValue=opts.Convert(min(im(:)));
    end
    if (isfield(info, 'LargestImagePixelValue'))
        info.LargestImagePixelValue=opts.Convert(max(im(:)));
    end
    if (isfield(info, 'RescaleIntercept'))
        info.RescaleIntercept=0;
    end
    if (isfield(info, 'RescaleSlope'))
        info.RescaleSlope=1;
    end
    info.WindowCenter=opts.Level;
    info.WindowWidth=opts.Window;
	info.SamplesPerPixel=1;
    info.BitsAllocated=16;
    info.BitsStored=16;
    info.HighBit=15;
    info.PixelRepresentation=flag_signed;
    info.SeriesDescription = opts.SeriesDescription;
    info.SeriesInstanceUID = opts.SeriesInstanceUID;
    info.SOPInstanceUID = dicomuid;
    if isempty(opts.SeriesNumber)
        opts.SeriesNumber=info.SeriesNumber*100;
    end
    info.SeriesNumber = opts.SeriesNumber;
    info.InstanceNumber = slice;
    outfile=fullfile(outdir,sprintf('IM%05d.dcm', slice));
    print_progress(outfile);
    dicomwrite(opts.Convert(M(:,:,slice))',outfile, ...
        'CreateMode', 'copy', 'WritePrivate', true, info);
end
fprintf('\n');


    function print_progress(arg)
        num=length(progress);
        num=num-numel(regexp(progress, '\\\\'));
        for ii=1:num; fprintf('\b'); end
        progress=['Writing file ' arg];
        progress=regexprep(progress,'\','\\\\');
        fprintf(progress);
    end

    function y=convert(x)
        if flag_signed
            y=int16(x*1000);
        else
            y=uint16(x*1000);
        end
    end
end