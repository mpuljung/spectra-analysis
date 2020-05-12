function pcf
%this command starts a gui to select *.spe files for analysis
%designate files for bleaching correction, as desired
%fits bleaching curve and prompts user to accept/reject fit
%if accepted, the remaining data are first analyzed according to spe2
%data are corrected
%final output is a *.xlsx file with the raw data, subtracted data, and
%bleaching corrected subtracted data.
close all force
%the new version doesn't automatically share variables from nested
%functions with the main function. Seems pretty fucked up  to me. Anyway,
%the work-around is to define the variables in the main function first as
%an empty set, so here goes.
image=[];
FP=[];
FPbkgd=[];
FPcorr=[];
x_axis=[];
exposure_time=[];
ROIcorr=[];
SUBheader=[];
RAWheader=[];
RAWdata=[];
PEAKheader=[];
PEAKdataname=[];
PEAKdata=[];
PEAKtopbottom=[];
BKGDtopbottom=[];
bleach_time_axis=[];
bleach_fit_data=[];
xData=[];
yData=[];
fitresult=[];

%%suppress warnings
warning('off','curvefit:prepareFittingData:sizeMismatch');
warning('off','MATLAB:xlswrite:AddSheet');
%%

[roi_filename,pathname]=uigetfile('.spe','Select Files for ROI',...
    'C:\Users\dpag0575\Documents\Data\Imaging\Pixis','MultiSelect','on');
cd(pathname)


%sort all files in folder
D = dir('*.spe');
% Engine
S = [D(:).datenum].'; % you may want to eliminate . and .. first.
[~,S] = sort(S);
S = {D(S).name}; % Cell array of names in order by datenum. 

%intersection of selected files and all files sorted by datenum
roi_sorted=intersect(S,roi_filename,'stable');

%enter name of output file (appended with xlsx)
xlsfile=input('Type name of output file>','s');
xlsfile=[xlsfile '.xlsx'];

%checks if file already exists and prompts you to give another name if it
%does
if exist([pathname xlsfile])>0
    xlsfile=input('That name is taken. Try another>','s');
    xlsfile=[xlsfile '.xlsx'];
else
end

[~,c]=size(xlsfile);
figname=xlsfile(1:c-5);
%%
%read in first two images, plot them, and allow for selection of background region 
readSPE(char(roi_sorted(1)))




%%

%readSPE.m
    function readSPE(dirPath,filename)
        %% Start of Code

% parse optional input
if nargin>1
    filename = strcat(dirPath,filename);
else
    filename = dirPath;
end

% Open the file
fd = fopen(filename,'r');
if(fd < 0)
    error('Could not open file, bad filename')
end

% Get the image dimensions:
stripDim = getData(fd, '2A', 'uint16');     %first dim
pixelDim = getData(fd, '290', 'uint16');    %second dim
nDim = getData(fd, '5A6', 'uint32');        %third dim

% Get the pixel data type
dataType = getData(fd, '6C', 'uint16');

% Get the image
fseek(fd, hex2dec('1004'), 'bof');



image = zeros([pixelDim,stripDim,nDim]);
switch dataType
    case 0     % single precision float (4 bytes)
        image = single(image);      %maintain datatype in function output
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'float32')';
        end
    
    case 1    % long int (4 bytes)
        image = int32(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int32')';
        end
    
    case 2    % short int (2 bytes)
        image = int16(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int16')';
        end
    
    case 3    % short unsigned int (2 bytes)
        image = uint16(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'uint16')';
        end

    otherwise
        image = -1;
end
    end
%%
%generate x axis
x=linspace(0,1399,1340);

%generate y axis
y=linspace(0,399,400);
%make xy grid
[X,Y]=meshgrid(x,y);

%make image data dp
intensity=double(image(:,:,1));

%make figure that plots bf and fluorescent images of pipetted side by side
%and allows for roi and background selection

roi_fig=figure('Name','Select ROI and Background','Position',[10 250 1500 500]);
%plot
subplot(1,2,1);
h=mesh(X(:,421:920),Y(:,421:920),intensity(:,421:920));
%change view
view(0,90)
title('Bright Field');
set(gca,'YDir','reverse')
pbaspect([500 400 1]);
axis tight
saveas(roi_fig,[figname 'ROI.png']);  

%read in first two images, plot them, and allow for selection of background region 

readSPE(char(roi_sorted(2)));
%make image data dp
intensity_fluor=double(image(:,:,1));

subplot(1,2,2);
h=mesh(X(:,421:920),Y(:,421:920),intensity_fluor(:,421:920));
%change view
view(0,90)
title('Fluorescence')
set(gca,'YDir','reverse')
pbaspect([500 400 1]);
axis tight

saveas(roi_fig,[figname 'ROI.png']);


%select top and bottom of ROI with cursor
datacursormode on
disp('Click top of the spectrum and press Enter.')
[~,top]=ginput;
roundtop=round(top);
refresh
disp('Click bottom of the spectrum and press Enter.')
[~,bottom]=ginput;
roundbottom=round(bottom);

refresh
disp('Click top of the background region and press Enter.')
[~,bkgdtop]=ginput;
if bkgdtop<1
    bkgdtop=1;
end
roundbkgdtop=round(bkgdtop);
refresh
disp('Click bottom of the background region and press Enter.')
[~,bkgdbottom]=ginput;
if bkgdbottom>1399
    bkgbottom=1399;
end
roundbkgdbottom=round(bkgdbottom);
refresh


close

%%
%find file for FP spectrum and plot
[FP_name,pathname]=uigetfile('.spe','Select File for FP Spectrum',pathname,'MultiSelect','on');
FPanal(char(FP_name))



%%
%FPanal.m
function FPanal(dirPath,filename)
%opens SPE file using readSPE.m at the beginning
%creates x values in units of wavelength, which are copied right out of the
%original spe file, as apparently the increments aren't linear.

%Data are read in as a matrix of intensity values.  Change to double
%precision.

%create 3D (mesh) plot and rotate it.



%% Start of Code


% parse optional input

if nargin>1
    filename = strcat(dirPath,filename);
else
    filename = dirPath;
end

% Open the file
fd = fopen(filename,'r');
if(fd < 0)
    error('Could not open file, bad filename')
end

% Get the image dimensions:
stripDim = getData(fd, '2A', 'uint16');     %first dim
pixelDim = getData(fd, '290', 'uint16');    %second dim
nDim = getData(fd, '5A6', 'uint32');        %third dim

% Get the pixel data type
dataType = getData(fd, '6C', 'uint16');

% Get the image
fseek(fd, hex2dec('1004'), 'bof');

image = zeros([pixelDim,stripDim,nDim]);
switch dataType
    case 0     % single precision float (4 bytes)
        image = single(image);      %maintain datatype in function output
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'float32')';
        end
    
    case 1    % long int (4 bytes)
        image = int32(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int32')';
        end
    
    case 2    % short int (2 bytes)
        image = int16(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int16')';
        end
    
    case 3    % short unsigned int (2 bytes)
        image = uint16(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'uint16')';
        end

    otherwise
        image = -1;
end

%BEGIN SAM USHER CODE
%offset of footer 
obj.FooterOffset = getData(fd, '2A6', 'uint64'); 

%navigating to start of the footer
fseek(fd, obj.FooterOffset, 'bof'); 

%reading footer in as characters
xml1 = fread(fd, Inf, '*char');
xml2 = char(xml1');

%splitting character into the xml segments
segments = strings(0);
remain = xml2;
while (remain ~= "")
[token,remain] = strtok(remain, '<>');
segments = [segments ; token];
end

%lengthy conversion of string to numeric vector.  This should definitely be
%shorter.
wavelengthraw = segments(10);

splitwavelengthraw = strsplit(wavelengthraw, ',');

needstransposing = str2double(splitwavelengthraw);

lambda = transpose(needstransposing);
%END SAM USHER CODE
%%
%generate y axis
y=linspace(0,399,400);
%make xy grid
[nm,Y]=meshgrid(lambda,y);

%make image data dp
intensity=double(image(:,:,1));

FP=mean(intensity(roundtop:roundbottom,400:1339));
FPbkgd=mean(intensity(roundbkgdtop:roundbkgdbottom,400:1339));
FPcorr=FP-FPbkgd;

x_axis=(lambda(400:1339));



%begin matrices/cell arrays for each excel sheet, titles, data
[~,c]=size(FP_name);
outputname=(FP_name(1:c-4));
FPheader=[{'nm'}, {[outputname 'RAW']},{[outputname 'BKGD']}];
FPdata=[x_axis,FP',FPbkgd'];


end

%%
%select files for analysis
[anal_name,pathname]=uigetfile('.spe','Select Files for Analysis',pathname,'MultiSelect','on');
[~,last]=size(anal_name);
cd(pathname)
%select files for bleaching correction
[bleach_name,pathname]=uigetfile('.spe','Select Files for Bleaching Correction',pathname,'MultiSelect','on');
    [~,bleach_sorted,~]=intersect(anal_name,bleach_name,'stable'); %intersection of selected bleaching files and all files selected for analysis
    bleach_intersect=intersect(bleach_name,anal_name,'stable'); %find the intersection between files selected for bleaching correction and all files being analyzed
    if isequal(bleach_name,bleach_intersect)==0 %see if the bleach intersect is the same as the files selected for bleaching
        display('Error: bleach correction files not included in analysis') %if they aren't, files were selected for bleaching correction outside of analysis files.
        return
    else 
    end
        
 


SPEanal(char(anal_name(1)))

    

%SPEanal.m
    function SPEanal(dirPath,filename)
%opens SPE file using readSPE.m at the beginning
%creates x values in units of wavelength, which are copied right out of the
%original spe file, as apparently the increments aren't linear.

%Data are read in as a matrix of intensity values.  Change to double
%precision.

%create 3D (mesh) plot and rotate it.



%% Start of Code


% parse optional input
if nargin>1
    filename = strcat(dirPath,filename);
else
    filename = dirPath;
end

% Open the file
fd = fopen(filename,'r');
if(fd < 0)
    error('Could not open file, bad filename')
end

% Get the image dimensions:
stripDim = getData(fd, '2A', 'uint16');     %first dim
pixelDim = getData(fd, '290', 'uint16');    %second dim
nDim = getData(fd, '5A6', 'uint32');        %third dim

% Get the pixel data type
dataType = getData(fd, '6C', 'uint16');

% Get the image
fseek(fd, hex2dec('1004'), 'bof');

image = zeros([pixelDim,stripDim,nDim]);
switch dataType
    case 0     % single precision float (4 bytes)
        image = single(image);      %maintain datatype in function output
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'float32')';
        end
    
    case 1    % long int (4 bytes)
        image = int32(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int32')';
        end
    
    case 2    % short int (2 bytes)
        image = int16(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int16')';
        end
    
    case 3    % short unsigned int (2 bytes)
        image = uint16(image);
        for i=1:nDim
            image(:,:,i) = fread(fd, [stripDim,pixelDim], 'uint16')';
        end

    otherwise
        image = -1;
end

%BEGIN SAM USHER CODE
%offset of footer 
obj.FooterOffset = getData(fd, '2A6', 'uint64'); 

%navigating to start of the footer
fseek(fd, obj.FooterOffset, 'bof'); 

%reading footer in as characters
xml1 = fread(fd, Inf, '*char');
xml2 = char(xml1');

%splitting character into the xml segments
segments = strings(0);
remain = xml2;
while (remain ~= "")
[token,remain] = strtok(remain, '<>');
segments = [segments ; token];
end

%lengthy conversion of string to numeric vector.  This should definitely be
%shorter.
wavelengthraw = segments(10);

splitwavelengthraw = strsplit(wavelengthraw, ',');

needstransposing = str2double(splitwavelengthraw);

lambda = transpose(needstransposing);
%END SAM USHER CODE
%%
%get exposure time from segments
[~,~,exposure_location]=intersect('ExposureTime type="Double"',segments);
exposure_time=str2num(segments(exposure_location+1))/1000;
%%
%make image data dp
intensity=double(image(:,:,1));

ROI=mean(intensity(roundtop:roundbottom,400:1339));
ROIbkgd=mean(intensity(roundbkgdtop:roundbkgdbottom,400:1339));
ROIcorr=ROI-ROIbkgd;

x_axis=(lambda(400:1339));



%find peak in data
win=31;
[~,length]=size(ROIcorr);
means=zeros((length-(win-1)),2);
means(:,1)=x_axis((((win-1)/2)+1):(length-((win-1)/2)));

x=((win-1)/2)+1;
y=1;

%calculate the mean of a window at each point 
%(minus the ends, i.e. first and last 15 points)
while x<=(length-((win-1)/2)) 
    means(y,2)=mean(ROIcorr((x-((win-1)/2)):(x+((win-1)/2))));
    x=x+1;
    y=y+1;
end
[max_val,center_pos]=max(means(:,2));   %gives you max and the index of the row where the max was found


%begin matrices/cell arrays for each excel sheet, titles, data
first_anal=anal_name(1);
[~,c]=size(filename);
outputname=(filename(1:c-4));
SUBheader=[{'nm'},{outputname}];
SUBdata=[x_axis,ROIcorr'];
RAWheader=[{'nm'}, {[outputname 'RAW']},{[outputname 'BKGD']}];
RAWdata=[x_axis,ROI',ROIbkgd'];
PEAKheader=[{' '},{'nm'},{'time'},{'peak'},{'corrected'},{'top_of_ROI'}, {'bottom_of_ROI'}];
PEAKdataname=[{outputname}];
PEAKdata=[means(center_pos,1),max_val];
PEAKtopbottom=[roundtop,roundbottom];
BKGDtopbottom=[roundbkgdtop,roundbkgdbottom];


[~,number_of_files]=size(anal_name);
for  inc=(2:1:number_of_files)
callnext(char(anal_name(inc)))
end

function callnext(dirPath,filename) %calls loads subsequent files, 
                    %applies the same ROIs and
                  %appends the same Excel sheet
    if nargin>1;
        filename = strcat(dirPath,filename);
    else
        filename = dirPath;
    end

    % Open the file
    fd = fopen(filename,'r');
    if(fd < 0);
        error('Could not open file, bad filename')
    end

    % Get the image dimensions:
    stripDim = getData(fd, '2A', 'uint16');     %first dim
    pixelDim = getData(fd, '290', 'uint16');    %second dim
    nDim = getData(fd, '5A6', 'uint32');        %third dim

    % Get the pixel data type
    dataType = getData(fd, '6C', 'uint16');

    % Get the image
    fseek(fd, hex2dec('1004'), 'bof');

    image = zeros([pixelDim,stripDim,nDim]);
    switch dataType
        case 0;     % single precision float (4 bytes)
            image = single(image);      %maintain datatype in function output
            for i=1:nDim;
                image(:,:,i) = fread(fd, [stripDim,pixelDim], 'float32')';
            end
    
        case 1;    % long int (4 bytes)
            image = int32(image);
            for i=1:nDim;
                image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int32')';
            end
    
        case 2;    % short int (2 bytes)
            image = int16(image);
            for i=1:nDim;
                image(:,:,i) = fread(fd, [stripDim,pixelDim], 'int16')';
            end
    
        case 3;    % short unsigned int (2 bytes)
            image = uint16(image);
            for i=1:nDim;
                image(:,:,i) = fread(fd, [stripDim,pixelDim], 'uint16')';
            end

        otherwise
            image = -1;
    end
    %make image data dp
    intensity=double(image);

    ROI=mean(intensity(roundtop:roundbottom,400:1339));
    ROIbkgd=mean(intensity(roundbkgdtop:roundbkgdbottom,400:1339));
    ROIcorr=ROI-ROIbkgd;
    newcenter=center_pos+((win-1)/2);
    startm=newcenter-((win-1)/2);
    finm=newcenter+((win-1)/2);
    newpeak=mean(ROIcorr(startm:finm));

    %append matrices/cell arrays for each excel sheet, titles, data
    [~,c]=size(filename);
    outputname=(filename(1:c-4));
    
    [~,c]=size(SUBheader);
    SUBheader(:,c+1)={outputname};
    [~,c]=size(SUBdata);
    SUBdata(:,c+1)=ROIcorr';
    
    [~,c]=size(RAWheader);
    RAWheader(:,c+1)={[outputname 'RAW']};
    RAWheader(:,c+2)={[outputname 'BKGD']};
    [~,c]=size(RAWdata);
    RAWdata(:,c+1)=ROI';
    RAWdata(:,c+2)=ROIbkgd';    
    
    [r,~]=size(PEAKdataname);
    PEAKdataname(r+1,:)=[{outputname}];
    [r,~]=size(PEAKdata);
    PEAKdata(r+1,:)=[means(center_pos,1),newpeak];   
    

end
    
%find data from PEAKdata corresponding to max of bleaching steps
bleach_time_axis=exposure_time*(bleach_sorted(:)-bleach_sorted(1))';
bleach_fit_data=[PEAKdata(bleach_sorted,2)/max(PEAKdata(bleach_sorted,2))];



%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'bleaching' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 14-Feb-2018 12:59:17


%% Fit: 'bleaching'.
[xData, yData] = prepareCurveData(bleach_time_axis, bleach_fit_data);

% Set up fittype and options.
ft = fittype( 'a*exp(-b*bleach_time_axis)+(1-a)', 'independent', 'bleach_time_axis', 'dependent', 'bleach_fit_data' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.655796065597017 0.69832204678041];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
bleaching_corr=coeffvalues(bleachingfit(bleach_time_axis,bleach_fit_data));

%time axis
[~,fileindex,~]=intersect(S,anal_name,'stable');
newtimeaxis=exposure_time*(fileindex(:)-fileindex(1)); %returns column vector with time points for every file in 'sorted' based on their position in S 

denom=bleaching_corr(1)*exp(-newtimeaxis*bleaching_corr(2))+(1-bleaching_corr(1));
corrPEAKdata=PEAKdata(:,2)./denom; %correct peak data for bleaching
[~,c]=size(SUBdata);
corrSUBdata=SUBdata(:,2:c)./denom';






accept=input('Accept fit (y/n)?','s');

while accept=='n'
    
    close
    %select new files for bleaching correction
    [bleach_name,pathname]=uigetfile('.spe','Select Files for Bleaching Correction',pathname,'MultiSelect','on');
    [~,bleach_sorted,~]=intersect(anal_name,bleach_name,'stable'); %intersection of selected bleaching files and all files selected for analysis
    bleach_intersect=intersect(bleach_name,anal_name,'stable'); %find the intersection between files selected for bleaching correction and all files being analyzed
    if isequal(bleach_name,bleach_intersect)==0 %see if the bleach intersect is the same as the files selected for bleaching
        display('Error: bleach correction files not included in analysis') %if they aren't, files were selected for bleaching correction outside of analysis files.
        return
    else 
    end
    %% Fit: 'bleaching'.
    
    bleach_time_axis=exposure_time*(bleach_sorted(:)-bleach_sorted(1))';
    bleach_fit_data=[PEAKdata(bleach_sorted,2)/max(PEAKdata(bleach_sorted,2))];
    [xData, yData] = prepareCurveData(bleach_time_axis, bleach_fit_data);

    % Set up fittype and options.
    ft = fittype( 'a*exp(-b*bleach_time_axis)+(1-a)', 'independent', 'bleach_time_axis', 'dependent', 'bleach_fit_data' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.655796065597017 0.69832204678041];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    % Plot fit with data.
bleaching_corr=coeffvalues(bleachingfit(bleach_time_axis,bleach_fit_data));

%time axis
[~,fileindex,~]=intersect(S,anal_name,'stable');
newtimeaxis=exposure_time*(fileindex(:)-fileindex(1)); %returns column vector with time points for every file in 'sorted' based on their position in S 

denom=bleaching_corr(1)*exp(-newtimeaxis*bleaching_corr(2))+(1-bleaching_corr(1));
corrPEAKdata=PEAKdata(:,2)./denom; %correct peak data for bleaching
[~,c]=size(SUBdata);
corrSUBdata=SUBdata(:,2:c)./denom';

accept=input('Accept fit (y/n)?','s');


end








    end

bleaching_corr=coeffvalues(bleachingfit(bleach_time_axis,bleach_fit_data));

%time axis
[~,fileindex,~]=intersect(S,anal_name,'stable');
newtimeaxis=exposure_time*(fileindex(:)-fileindex(1)); %returns column vector with time points for every file in 'sorted' based on their position in S 

denom=bleaching_corr(1)*exp(-newtimeaxis*bleaching_corr(2))+(1-bleaching_corr(1));
corrPEAKdata=PEAKdata(:,2)./denom; %correct peak data for bleaching
[~,c]=size(SUBdata);
corrSUBdata=SUBdata(:,2:c)./denom';

outputfig=figure('Name',figname,'Position',[500 80 800 700]);
subplot(3,2,1)
h = plot(fitresult, xData, yData,'o');
ylim([0 1])
legend( h, 'data', 'fit', 'Location', 'SouthEast' );
% Label axes
title('bleaching')
xlabel 'time s'
ylabel 'normalized fluorescence'
grid on

% Plot corrected bleaching data

subplot(3,2,2)
bleachplot = plot(x_axis,corrSUBdata(:,bleach_sorted));
% Label axes
title('bleaching corrected')
xlabel 'nm'
ylabel 'normalized fluorescence'
axis tight



CORRECTEDheader1=[{'y=a*(-x*k)+(1-a)'}];
CORRECTEDheader2=[{'time s'},{'raw'},{'norm'},{'corr'},{'amp'},{'1/tau'}];


%find names and indices of files that were analyzed that were NOT part of
%bleaching correction
[data_used,data_index]=setdiff(anal_name,bleach_intersect','stable');
%Now assemble all the bits and write a new xlsx

xlswrite(xlsfile,[{'nm'},{'FP fluor'},{'FP bkgd'},{'FP sub'}],'FP','A1');
xlswrite(xlsfile,x_axis,'FP','A2');
xlswrite(xlsfile,FP','FP','B2');
xlswrite(xlsfile,FPbkgd','FP','C2');
xlswrite(xlsfile,FPcorr','FP','D2');
xlswrite(xlsfile,SUBheader,'Subtracted','A1');
xlswrite(xlsfile,SUBdata,'Subtracted','A2');
xlswrite(xlsfile,RAWheader,'Raw','A1');
xlswrite(xlsfile,RAWdata,'Raw','A2');
xlswrite(xlsfile,PEAKheader,'Peaks','A1');
xlswrite(xlsfile,PEAKtopbottom,'Peaks','F2');
xlswrite(xlsfile,BKGDtopbottom,'Peaks','F3');
xlswrite(xlsfile,PEAKdataname,'Peaks','A2');
xlswrite(xlsfile,PEAKdata(:,1),'Peaks','B2');
xlswrite(xlsfile,newtimeaxis,'Peaks','C2');
xlswrite(xlsfile,PEAKdata(:,2),'Peaks','D2');
xlswrite(xlsfile,corrPEAKdata,'Peaks','E2');
xlswrite(xlsfile,CORRECTEDheader1,'bleaching','A1');
xlswrite(xlsfile,CORRECTEDheader2,'bleaching','A2');
xlswrite(xlsfile,bleach_time_axis','bleaching','A3');
xlswrite(xlsfile,PEAKdata(bleach_sorted,2),'bleaching','B3');
xlswrite(xlsfile,PEAKdata(bleach_sorted,2)/max(PEAKdata(bleach_sorted,2)),'bleaching','C3');
xlswrite(xlsfile,corrPEAKdata(bleach_sorted),'bleaching','D3');
xlswrite(xlsfile,bleaching_corr,'bleaching','E3');
xlswrite(xlsfile,denom','corrected','B1');
xlswrite(xlsfile,SUBheader,'corrected','A2');
xlswrite(xlsfile,x_axis,'corrected','A3');
xlswrite(xlsfile,corrSUBdata,'corrected','B3');
xlswrite(xlsfile,{'file' 'F/Fmax'},'Peaks','I1');
xlswrite(xlsfile,[anal_name(data_index(1)-1); data_used'],'Peaks','I2');
xlswrite(xlsfile,[corrPEAKdata((data_index(1)-1))/corrPEAKdata((data_index(1)-1)); corrPEAKdata(data_index)./corrPEAKdata((data_index(1)-1))],'Peaks','J2');


% Plot FP spectrum and first ANAP spectrum
subplot(3,2,3)
FPplot = plot(x_axis,[FPcorr' SUBdata(:,2)]);
% Label axes
title('FP spectrum')
xlabel 'nm'
ylabel 'fluorescence'
legend(FPplot,'FP', 'ANAP', 'Location', 'NorthEast' );
axis tight


subplot(3,2,4)
concatdata=[corrSUBdata(:,data_index(1)-1) corrSUBdata(:,data_index)];
specplot = plot(x_axis,concatdata);
% Label axes
title('corrected spectra')
xlabel 'nm'
ylabel 'fluorescence'
leg=[anal_name(data_index(1)-1) data_used];
legend(specplot, leg, 'Location', 'NorthEast' )
axis tight


subplot(3,2,[5,6])
x1=categorical([anal_name(data_index(1)-1); data_used']);
x1=reordercats(x1,[anal_name(data_index(1)-1); data_used']);
y1=[corrPEAKdata((data_index(1)-1))/corrPEAKdata((data_index(1)-1)); corrPEAKdata(data_index)./corrPEAKdata((data_index(1)-1))];
summaryplot=plot(x1,y1,'o');
% Label axes
ylabel 'F/Fmax'
ylim([-0.1 1.1])


saveas(outputfig,[figname '.png']);


beep
  
end

%% 
% getData() reads one piece of data at a specific location
% 
function data = getData(fd, hexLoc, dataType)
% Inputs: fd - int    - file descriptor
%     hexLoc - string - location of data relative to beginning of file
%   dataType - string - type of data to be read
%
fseek(fd, hex2dec(hexLoc), 'bof');
data = fread(fd, 1, dataType);

end
