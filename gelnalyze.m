function gelnalyze(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A=imread(filename); %read in file
[r,c]=size(A);
grid=ones(r,c); %make a grid for mesh plot
figure
mesh(grid,A)%plot
view(0,-90) %change view so it's the same as originial 2D image

%start with lane 1
n=0;
lane=num2str(n+1);

%pick ROI
datacursormode on
display('click left side of lane and press enter')
[left,~]=ginput;
left=round(left);
refresh
display('click right side of lane and press enter')
[right,~]=ginput;
right=round(right);

%calculate mean across columns for all rows
ROI=A(:,left:right);
ROI_out=mean(ROI,2);

%pick bkgd
refresh
display('click left side of background and press enter')
[left_bkgd,~]=ginput;
left_bkgd=round(left_bkgd);
refresh
display('click right side of background and press enter')
[right_bkgd,~]=ginput;
right_bkgd=round(right_bkgd);

%calculate bkgd across columns for all rows
BKGD=A(:,left_bkgd:right_bkgd);
BKGD_out=mean(BKGD,2);

%subtract bkgd from mean
SUB_out=ROI_out-BKGD_out;

%start cell array for header
header(1,3*n+1)={['ROI lane' lane]};
header(1,3*n+2)={['BKGD lane' lane]};
header(1,3*n+3)={['SUB lane' lane]};

%start array for data
data(1:r,3*n+1)=ROI_out;
data(1:r,3*n+2)=BKGD_out;
data(1:r,3*n+3)=SUB_out;

prompt='pick more lanes? (y/n)';
str=input(prompt, 's');

while str=='y'
    
n=n+1;
lane=num2str(n+1);

%pick nextROI
datacursormode on
display('click left side of lane and press enter')
[left,~]=ginput;
left=round(left);
refresh
display('click right side of lane and press enter')
[right,~]=ginput;
right=round(right);

%calculate mean across columns for all rows
ROI=A(:,left:right);
ROI_out=mean(ROI,2);

%pick bkgd
refresh
display('click left side of background and press enter')
[left_bkgd,~]=ginput;
left_bkgd=round(left_bkgd);
refresh
display('click right side of background and press enter')
[right_bkgd,~]=ginput;
right_bkgd=round(right_bkgd);

%calculate bkgd across columns for all rows
BKGD=A(:,left_bkgd:right_bkgd);
BKGD_out=mean(BKGD,2);

%subtract bkgd from mean
SUB_out=ROI_out-BKGD_out;

%add to cell array for header
header(1,3*n+1)={['ROI lane' lane]};
header(1,3*n+2)={['BKGD lane' lane]};
header(1,3*n+3)={['SUB lane' lane]};

%add to array for data
data(1:r,3*n+1)=ROI_out;
data(1:r,3*n+2)=BKGD_out;
data(1:r,3*n+3)=SUB_out;

prompt='pick more lanes? (y/n)';
str=input(prompt, 's');

end


%prepare *.xlsx file name
[~,letters]=size(filename);
xlsfile=[filename(1:letters-3) 'xlsx'];


%write *.xlsx file
%headers

xlswrite(xlsfile,header,'Sheet1','A1');
xlswrite(xlsfile,data,'Sheet1','A2');

beep
end

