function spegui
%this command starts a gui to select *.spe files for analysis
%designate files for bleaching correction, as desired
%fits bleaching curve and prompts user to accept/reject fit
%if accepted, the remaining data are first analyzed according to spe2
%data are corrected
%final output is a *.xlsx file with the raw data, subtracted data, and
%bleaching corrected subtracted data.
close all force
[filename,pathname]=uigetfile('.spe','Select Files to Analyze','C:\Users\dpag0575\Documents\Data\Imaging\Pixis','MultiSelect','on');
cd(pathname);

%sort all files in folder
D = dir('*.spe');
% Engine
S = [D(:).datenum].'; % you may want to eliminate . and .. first.
[~,S] = sort(S);
S = {D(S).name}; % Cell array of names in order by datenum. 

%intersection of selected files and all files sorted by datenum
sorted=intersect(S,filename,'stable');
[~,number_of_files]=size(sorted);

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

%The next bit puts up a selection box to pick bleach correction files.
%Files may be picked in any order, as it will sort them into chronological
%order. A button must be pressed to select.
%The final cell array with file names is bleach_corr_sorted
function getbleach
        bleach_corr=get(lbx,'Value');
        bleach_corr_sorted=intersect(sorted,bleach_corr,'stable');
    end
newfig=uifigure('Name','Select Bleaching Correction Files','Position',[50 50 500 275]);
lbx=uilistbox(newfig,'Position',[25 62.5 450 200],'Items',sorted,'Multiselect','on');
pb=uibutton(newfig,'push','Text','select','Position',[225 20 50 20],'ButtonPushedFcn', @(btn,event) getbleach);
    


SPEanal(char(sorted(1)))

%SPEanal.m
    function SPEanal(dirPath,filename)
%opens SPE file using readSPE.m at the beginning
%creates x values in units of wavelength, which are copied right out of the
%original spe file, as apparently the increments aren't linear.

%Data are read in as a matrix of intensity values.  Change to double
%precision.

%create 3D (mesh) plot and rotate it.

% readSPE.m  Read princeton instrument *.SPE image file
%
% image = SPEanal(filename)
%       = SPEanal(dirPath, filename)
%
% where: image - 3D array of image(s)
%     filename - path and filename (string)
%      dirPath - optional directory path (string)
% 
% Image is returned as a 3D array, where each image is the first 2
% dimensions, and successive images are stored along the 3rd dimension.
%
% The image is stored as it is shown in WinVIEW; the first two dimensions 
% are stored as [pixel,stripe]
% 
% returns -1 if unsuccessful, though if any of the low level access
% functions error out, this method will just fail with their error
% descriptions. Only the uint16 pixel datatype feature has been tested, so
% be more cautious than usual if using other datatypes.
%
% There is far more information stored in the SPE header (4100 bytes).
% This simple method only pulls the information that is required for the 
% author's data processing. For further documentation of the header byte 
% structure, seek the sources I used:
%  
%  Matt Clarke's Python Script Description:
%   http://www.scipy.org/Cookbook/Reading_SPE_files
%
%  Stuart B. Wilkins's Python Code:
%   https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/ccd/files.py 
%
%  ImageJ's plugin to import SPE files, written in Java:
%   http://rsbweb.nih.gov/ij/
%  

% Author: Carl Hall
% Date: June 2012

% Modification History
%  April 2012 - Original Code
%   June 2012 - Modified to maintain output array as same datatype as SPE
%               file

%% SPE Header byte Structure (not fully checked)
% 
% Offset    Bytes    Type     Description
% 0x0000    2        int      Controller Version
% 0x0002    2        int      Logic Output
% 0x0004    2        int      AppHiCapLowNoise
% 0x0006    2        int      dxDim
% 0x0008    2        int      timingMode
% 0x000A    4        float    Exposure
% 0x000E    2        int      vxDim
% 0x0010    2        int      vyDim
% 0x0012    2        int      dyDim
% 0x0014    10       string   Date
% 0x0024    4        float    DetTemperature
% 0x0028    2        int      DectectorType
% 0x002A    2        int      xDim
% 0x002C    2        int      TriggerDiode
% 0x002E    4        float    DelayTime
% 0x0032    2        int      ShutterControl
% 0x0034    2        int      AbsorbLive
% 0x0036    2        int      AbsorbMode
% 0x0038    2        int      canDoVirtualChip
% 0x003A    2        int      thresholdMinLive
% 0x003C    4        float    threshholdMin
% 0x0040    2        int      thresholdMaxLive
% 0x0042    4        float    threshholdMax
% 0x006C    2        int      Data Type
% 0x00AC    7        string   Time
% 0x00BC    2        int      ADCOffset
% 0x00BE    2        int      ADCRate
% 0x00C0    2        int      ADCType
% 0x00C2    2        int      ADCRes
% 0x00C4    2        int      ADCBitAdj
% 0x00C6    2        int      Gain
% 0x00C8    80       string   Comments
% 0x0118    80       string   Comments
% 0x0168    80       string   Comments
% 0x01B8    80       string   Comments
% 0x0208    80       string   Comments
% 0x0258    2        int      GeometricOps
% 0x0290    2        int      ydim
% 0x05A6    2        uint32   zdim
% 0x05D0    2        int      NumROIExperiment
% 0x05D2    2        int      NumROI
% 0x05D4    60       int[]    allROI    
% 
% 0x1004    ...      datatype   Image Data

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

%generate y axis
y=linspace(0,399,400);
%make xy grid
[nm,Y]=meshgrid(lambda,y);

%make image data dp
intensity=double(image(:,:,1));

%plot
figure
mesh(nm(:,400:1339),Y(:,400:1339),intensity(:,400:1339))
%change view
view(0,90)
xlabel('nm');
title(filename);
set(gca,'YDir','reverse')
axis tight

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

ROI=mean(intensity(roundtop:roundbottom,400:1339));
ROIbkgd=mean(intensity(roundbkgdtop:roundbkgdbottom,400:1339));
ROIcorr=ROI-ROIbkgd;

x_axis=(lambda(400:1339));

%find peaks in data
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

for  inc=(2:1:number_of_files)
callnext(char(sorted(inc)))
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
    end

    
      
%find data from PEAKdata corresponding to max of bleaching steps
[~,bleach_index,~]=intersect(sorted,bleach_corr_sorted,'stable'); %returns the location (number of elements in sorted) and therfore PEAKdata corresponding to bleach corr files
bleach_time_axis=10*(bleach_index(:)-bleach_index(1))';
bleach_fit_data=[PEAKdata(bleach_index,2)/max(PEAKdata(bleach_index,2))];


function [fitresult, gof] = bleachingfit(x, y)
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
[~,c]=size(xlsfile);
figname=xlsfile(1:c-5);
outputfig=figure('Name',figname,'Position',[500 80 800 700]);
subplot(2,2,1)
h = plot( fitresult, xData, yData,'o');
legend( h, 'data', 'fit', 'Location', 'NorthEast' );
% Label axes
title('bleaching')
xlabel 'time s'
ylabel 'normalized fluorescence'
grid on

end
bleaching_corr=coeffvalues(bleachingfit(bleach_time_axis,bleach_fit_data));


%time axis
[~,fileindex,~]=intersect(S,sorted,'stable');
newtimeaxis=10*(fileindex(:)-fileindex(1)); %returns column vector with time points for every file in 'sorted' based on their position in S 

denom=bleaching_corr(1)*exp(-newtimeaxis*bleaching_corr(2))+(1-bleaching_corr(1));
corrPEAKdata=PEAKdata(:,2)./denom; %correct peak data for bleaching
[~,c]=size(SUBdata);
corrSUBdata=SUBdata(:,2:c)./denom';

% Plot corrected bleaching data
[r,~]=size(bleach_fit_data);
subplot(2,2,2)
bleachplot = plot(x_axis,corrSUBdata(:,2:(r)));
% Label axes
title('bleaching corrected')
xlabel 'nm'
ylabel 'normalized fluorescence'
axis tight

concaxis=[5e-8;5e-7;5e-6;5e-5;5e-4;1e-3];
[rr,~]=size(corrPEAKdata);
concdata=1-corrPEAKdata(r+2:2:rr)/corrPEAKdata(r+1);



CORRECTEDheader1=[{'y=a*(-x*k)+(1-a)'}];
CORRECTEDheader2=[{'time s'},{'raw'},{'norm'},{'corr'},{'amp'},{'1/tau'}];
[k,~]=size(bleach_time_axis);

%Now assemble all the bits and write a new xlsx
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
xlswrite(xlsfile,PEAKdata(bleach_index,2),'bleaching','B3');
xlswrite(xlsfile,PEAKdata(bleach_index,2)/max(PEAKdata(bleach_index,2)),'bleaching','C3');
xlswrite(xlsfile,corrPEAKdata(1:k),'bleaching','D3');
xlswrite(xlsfile,bleaching_corr,'bleaching','E3');
xlswrite(xlsfile,denom','corrected','B1');
xlswrite(xlsfile,SUBheader,'corrected','A2');
xlswrite(xlsfile,x_axis,'corrected','A3');
xlswrite(xlsfile,corrSUBdata,'corrected','B3');
xlswrite(xlsfile,{'conc (M)' '1-F/Fmax'},'Peaks','I1');
xlswrite(xlsfile,concaxis,'Peaks','I2');
xlswrite(xlsfile,concdata,'Peaks','J2');

plotspecdata=cat(2,corrSUBdata(:,r+1),corrSUBdata(:,r+2:2:rr));
subplot(2,2,3)
specplot = plot(x_axis,plotspecdata);
% Label axes
title('corrected spectra')
xlabel 'nm'
ylabel 'normalized fluorescence'
leg=[0;concaxis];
legend(specplot, num2str(leg), 'Location', 'NorthEast' )
axis tight
hold

subplot(2,2,4)
concrespplot=semilogx(concaxis,concdata,'o');
% Label axes
title('concentration-response')
xlabel '[nuc] (M)'
ylabel '1-F/Fmax'
axis tight
hold

saveas(outputfig,[figname '.png']);
saveas(outputfig,[figname '.pdf']);

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

