function [idata,info]=ParRecReaderR5(PathName,FileName,Complex)

parrawdata = textscan(fopen([PathName,'\',FileName,'.PAR']),'%s', 'delimiter', '\n');
pardata=parrawdata{1};

for cntr=1:1000
    temp = strread(pardata{cntr},'%s', 'delimiter', ':');
    
    %**********************************************************************
    % Detection of the end of the header
    %**********************************************************************
    if isempty(temp)
        break;
    end
    
    temp2 = temp{1};
    if (temp2(1) == '#')
        continue % Ignore comments
    end    
    %**********************************************************************
    % Number of slices
    %**********************************************************************
    if (strfind(temp2,'Max. number of slices/locations')>1)
        info.NSlices = str2num(temp{2});
    end
    %**********************************************************************
    % Inplane resolution (Acquisition)
    %**********************************************************************
    if (strfind(temp2,'Scan resolution  (x, y)')>1)
        [Nx Ny] = strtok(temp{2},'  ');
        info.Nx=str2num(Nx);
        info.Ny=str2num(Ny);
    end
    %**********************************************************************
    % Echos
    %**********************************************************************
    if (strfind(temp2,'Max. number of echoes')>1)
        info.NEchos = str2num(temp{2});
    end
    %**********************************************************************
    % Cardiac phases
    %**********************************************************************
    if (strfind(temp2,'Max. number of cardiac phases')>1)
        info.NCard = str2num(temp{2});
    end
    %**********************************************************************
    % Dynamics
    %**********************************************************************
    if (strfind(temp2,'Max. number of dynamics')>1)
        info.NDyns = str2num(temp{2});
    end
    %**********************************************************************
    % Mixes
    %**********************************************************************
    if (strfind(temp2,'Max. number of mixes')>1)
        info.NMixes = str2num(temp{2});
    end
    
      %Angle (en radians)
    if (strfind(temp2,'Angulation midslice(ap,fh,rl)[degr]')>1)
        [rotY remains] = strtok(temp{2},'  ');
        [rotZ rotX]  = strtok(remains);
        info.rotY   =str2num(rotY)*pi/180;
        info.rotZ   =str2num(rotZ)*pi/180;
        info.rotX   =str2num(rotX)*pi/180;
    end
  
end
cntr=cntr+1;

temp   = strread(pardata{cntr},'%s');

%**************************************************************************
% Geometrical information
%**************************************************************************
info.xRes   = str2num(temp{10});
info.yRes   = str2num(temp{11});
info.zThick = str2num(temp{23});
info.zSpace = str2num(temp{24});
info.xSpace = str2num(temp{29});
info.ySpace = str2num(temp{30});

%**************************************************************************
% Magnitude scaling factors
% %**************************************************************************
RI(1)  = str2num(temp{12});
RS(1)  = str2num(temp{13});
SS(1)  = str2num(temp{14});

%**************************************************************************
% Phase scaling factors
%**************************************************************************
if (Complex==1)
    
    temp   = strread(pardata{cntr+1},'%s');
    RI(2)  = str2num(temp{12});
    RS(2)  = str2num(temp{13});
    SS(2)  = str2num(temp{14});
    
end
%**************************************************************************
% Read image
%**************************************************************************
recrawdata      = fread(fopen([PathName,'\',FileName,'.REC'],'r'),'uint16');

%**************************************************************************
% We make sure to take into account the different echos
%**************************************************************************
if (Complex==1)
    idata         = reshape(recrawdata,info.xRes,info.yRes,2,info.NEchos,info.NSlices); 
    idata=((idata(:,:,1,:,:)*RS(1)+RI(1))./(RS(1).*SS(1))).*exp(1i*(idata(:,:,2,:,:)*RS(2)+RI(2))./(RS(2).*SS(2)));
    idata=permute(squeeze(idata),[1 2 4 3]);
    idata_temp=zeros(size(idata));
    idata_temp(:,:,1:end/2,:)=idata(:,:,end/2+1:end,:);
    idata_temp(:,:,end/2+1:end,:)=idata(:,:,1:end/2,:);
    idata=idata_temp;
    clear idata_temp;
    info.TE=zeros(info.NEchos,1);
    for jj=1:info.NEchos;
        %**********************************************************************
        % We read the echo times
        %**********************************************************************
        temp   = strread(pardata{cntr},'%s');
        info.TE(jj)=str2num(temp{31})*0.001; % In second
        cntr=cntr+2;
          
    end   
else
    idata         = (reshape(recrawdata,info.xRes,info.yRes,info.NEchos,info.NSlices).*RS(1)+RI(1))./(RS(1).*SS(1));
    idata=permute(idata,[1 2 4 3]);
    idata_temp=zeros(size(idata));
    idata_temp(:,:,1:end/2,:)=idata(:,:,end/2+1:end,:);
    idata_temp(:,:,end/2+1:end,:)=idata(:,:,1:end/2,:);
    idata=idata_temp;
    clear idata_temp;
    info.TE=zeros(info.NEchos,1);
    for jj=1:info.NEchos;
        %**********************************************************************
        % We read the echo times
        %**********************************************************************
        temp   = strread(pardata{cntr},'%s');
        info.TE(jj)=str2num(temp{31})*0.001; % In second
        cntr=cntr+1;
          
    end
    
end
end
