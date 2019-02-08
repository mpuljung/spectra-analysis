%Copyright 2019 Michael Puljung, Samuel Usher

   %Licensed under the Apache License, Version 2.0 (the "License");
   %you may not use this file except in compliance with the License.
   %You may obtain a copy of the License at

    % http://www.apache.org/licenses/LICENSE-2.0

   %Unless required by applicable law or agreed to in writing, software
   %distributed under the License is distributed on an "AS IS" BASIS,
   %WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   %See the License for the specific language governing permissions and
   %limitations under the License.
   
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
xlsfile=[xlsfile ' washout.xlsx'];
time_inc=input('Enter time increment (in s)>');

%checks if file already exists and prompts you to give another name if it
%does
if exist([pathname xlsfile])>0
    xlsfile=input('That name is taken. Try another>','s');
    xlsfile=[xlsfile '.xlsx'];
else
end

 


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

lambda = str2double(splitwavelengthraw)';

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
height=roundbottom-roundtop;
refresh

ROI=mean(intensity(roundtop:roundbottom,400:1339));
ROIbkgd=mean(intensity(1:(height+1),400:1339));
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
ANAPlambda=find(floor(x_axis)>469); %finds location of all values in lambda that 
                                %round down to 470 or greater
ANAPnm=ANAPlambda(1);             %the first value
ANAPpeak=mean(ROIcorr((ANAPnm-15):(ANAPnm+15)));

%begin matrices/cell arrays for each excel sheet, titles, data
[~,c]=size(filename);
outputname=(filename(1:c-4));
SUBheader=[{'nm'},{outputname}];
SUBdata=[x_axis,ROIcorr'];
RAWheader=[{'nm'}, {[outputname 'RAW']},{[outputname 'BKGD']}];
RAWdata=[x_axis,ROI',ROIbkgd'];
PEAKheader=[{' '},{'nm'},{'time'},{'TNP'},{'ANAP'},{'subtracted&norm'},{'bkgd'},{'top_of_ROI'}, {'bottom_of_ROI'},{'peak/530'}];
PEAKdataname=[{outputname}];
PEAKdata=[means(center_pos,1),max_val];
ANAPdata=ANAPpeak;
PEAKtopbottom=[roundtop,roundbottom];


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
    ROIbkgd=mean(intensity(1:(height+1),400:1339));
    ROIcorr=ROI-ROIbkgd;
    newcenter=center_pos+((win-1)/2);
    startm=newcenter-((win-1)/2);
    finm=newcenter+((win-1)/2);
    newpeak=mean(ROIcorr(startm:finm));
    newANAPpeak=mean(ROIcorr((ANAPnm-15):(ANAPnm+15)));

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
    ANAPdata(r+1,:)=newANAPpeak;

    
    

end
    end
%time axis
[~,fileindex,~]=intersect(S,sorted,'stable');
newtimeaxis=time_inc*(fileindex(:)-fileindex(1)); %returns column vector with time points for every file in 'sorted' based on their position in S 

%choose time points after decay has gone 95% of the way for fitting a
%linear equation and correcting the washout data
[r,~]=size(PEAKdata);
mag=PEAKdata(1,2)-.99*(PEAKdata(1,2)-PEAKdata(r,2));
fitpoints=find(PEAKdata(:,2)<=mag);
testing=PEAKdata(fitpoints,2);
fit=polyfit(newtimeaxis(fitpoints),PEAKdata(fitpoints,2),1);
CORRdata=PEAKdata(:,2)-(fit(1)*newtimeaxis+fit(2));
NORMdata=CORRdata/CORRdata(1);

 % Plots
[~,c]=size(xlsfile);
figname=xlsfile(1:c-5);
outputfig=figure('Name',[figname 'wash'],'Position',[500 80 800 700]);  

      


% Plot subtracted data
subplot(2,2,1)
[~,c]=size(SUBdata);
corrplot = plot(x_axis,SUBdata(:,2:c));
% Label axes
title('subtracted spectra')
xlabel 'nm'
ylabel 'fluorescence'
axis tight
hold
[rend,~]=size(PEAKdata);
plot(means(center_pos,1),PEAKdata(1:rend,2),'mo')
hold

% Plot only first and last
subplot(2,2,2)
[~,c]=size(SUBdata);
firstlast(:,1)=SUBdata(:,2);
firstlast(:,2)=SUBdata(:,c);
firstplot = plot(x_axis,firstlast);
firstleg=[{'First'},{'Last'}];
legend(firstplot, firstleg, 'Location', 'NorthEast' );
% Label axes
title('subtracted spectra')
xlabel 'nm'
ylabel 'fluorescence'
axis tight
hold

%plot time course
subplot(2,2,3)
TIMEdata(:,1)=PEAKdata(:,2);
TIMEdata(:,2)=ANAPdata;
timeplot=plot(newtimeaxis,TIMEdata,'o');
% Label axes
title('TNP dissociation')
xlabel 'time (s)'
ylabel 'F/Fmax'
axis auto
hold
timeleg={'TNP peak' 'ANAP peak'};
legend(timeplot, timeleg, 'Location', 'NorthEast' );
axis auto
fity=fit(1)*newtimeaxis+fit(2);
plot(newtimeaxis(fitpoints),fity(fitpoints),'color',[1 0 0]);

%plot corrected time course
subplot(2,2,4)
corr_time=plot(newtimeaxis,CORRdata,'o');
title('corrected')
xlabel 'time (s)'
ylabel 'corrected F/Fmax'
axis auto
hold

outputfig_bkgd=figure('Name',[figname 'wash background'],'Position',[500 80 800 350]); 

%plot background spectra
subplot(1,2,1)
[~,c]=size(RAWdata);
corrplot = plot(x_axis,RAWdata(:,3:2:c));
% Label axes
title('background spectra')
xlabel 'nm'
ylabel 'fluorescence'
axis tight
hold

%plot background washout
subplot(1,2,2);
bkgdpeak=mean(RAWdata(353:383,3:2:c));
bkgd_time=plot(newtimeaxis, bkgdpeak,'o');
title('background wash')
xlabel 'time (s)'
ylabel 'F530/Fmax530'
axis auto
hold

%calculate peak/530 nm for last image
[r,~]=size(ANAPdata);
lastpeak=ANAPdata(r);
[~,c]=size(SUBdata);
lastTNPpeak=mean(SUBdata(353:383,c));
ratio=lastTNPpeak/lastpeak;



%Now assemble all the bits and write a new xlsx
xlswrite(xlsfile,SUBheader,'Subtracted','A1');
xlswrite(xlsfile,SUBdata,'Subtracted','A2');
xlswrite(xlsfile,RAWheader,'Raw','A1');
xlswrite(xlsfile,RAWdata,'Raw','A2');
xlswrite(xlsfile,PEAKheader,'Peaks','A1');
xlswrite(xlsfile,bkgdpeak','Peaks','G2');
xlswrite(xlsfile,PEAKtopbottom,'Peaks','H2');
xlswrite(xlsfile,PEAKdataname,'Peaks','A2');
xlswrite(xlsfile,PEAKdata(:,1),'Peaks','B2');
xlswrite(xlsfile,newtimeaxis,'Peaks','C2');
xlswrite(xlsfile,PEAKdata(:,2),'Peaks','D2');
xlswrite(xlsfile,ANAPdata,'Peaks','E2');
xlswrite(xlsfile,NORMdata,'Peaks','F2');
xlswrite(xlsfile,ratio,'Peaks','J2');



saveas(outputfig,[figname '.png']);
saveas(outputfig,[figname '.pdf']);
saveas(outputfig_bkgd,[figname 'bkgd.png']);
saveas(outputfig_bkgd,[figname 'bkgd.pdf']);

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

