function [idata,info]=ParRecReader(PathName,FileName,Complex)

% parrawdata = textscan(fopen([PathName,'\',FileName,'.PAR']),'%s', 'delimiter', '\n');
parrawdata = textscan(fopen([PathName,FileName,'.PAR']),'%s', 'delimiter', '\n');
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
    
    temp   = strread(pardata{cntr+info.NSlices*info.NEchos},'%s');
    RI(2)  = str2num(temp{12});
    RS(2)  = str2num(temp{13});
    SS(2)  = str2num(temp{14});
    
end
%**************************************************************************
% Read image
%**************************************************************************
% recrawdata      = fread(fopen([PathName,'\',FileName,'.REC'],'r'),'uint16');
recrawdata      = fread(fopen([PathName,FileName,'.REC'],'r'),'uint16');
%**************************************************************************
% We make sure to take into account the different echos
%**************************************************************************
if (Complex==1)
    idata_temp          = reshape(recrawdata,info.xRes,info.yRes,info.NSlices,info.NEchos*2);
    idata               = zeros(info.xRes,info.yRes,info.NSlices,info.NEchos);
    
    info.TE=zeros(info.NEchos,1);
    for jj=1:info.NEchos;
        %**********************************************************************
        % We read the echo times
        %**********************************************************************
        temp   = strread(pardata{cntr},'%s');
        info.TE(jj)=str2num(temp{31})*0.001; % In second
        cntr=cntr+info.NSlices;
        %**********************************************************************
        % We store the complex data
        %**********************************************************************
        idata(:,:,:,jj) = ( (idata_temp(:,:,:,jj).*RS(1)+RI(1))./(RS(1).*SS(1)) ).* ...
            (   cos( (idata_temp(:,:,:,jj+info.NEchos).*RS(2)+RI(2))./(RS(2).*SS(2)) )  + 1i.* ...
            sin( (idata_temp(:,:,:,jj+info.NEchos).*RS(2)+RI(2))./(RS(2).*SS(2)) ) );
        
    end
else
    idata         = (reshape(recrawdata,info.xRes,info.yRes,info.NSlices,info.NEchos).*RS(1)+RI(1))./(RS(1).*SS(1));
   
    info.TE=zeros(info.NEchos,1);
    for jj=1:info.NEchos;
        %**********************************************************************
        % We read the echo times
        %**********************************************************************
        temp   = strread(pardata{cntr},'%s');
        info.TE(jj)=str2num(temp{31})*0.001; % In second
        cntr=cntr+info.NSlices;
          
    end
    
end
end
