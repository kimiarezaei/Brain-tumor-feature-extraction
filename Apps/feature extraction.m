clear all
close all
clc
% preprocessing the images and extracting the features

% Read dicom format images
Files = dir('YOURFILEPATH\*.dcm');  % the folder in which your images exists
t = zeros (length(Files),7);

for i = 1 : length(Files)
filename = strcat('YOURFILEPATH\',Files(i).name);
image{i}=dicomread(filename);
%removing noise
image{i} = wiener2(image{i},[3 3]);
image{i} = medfilt2(image{i},[3 3]);
%descrete wavelet decomposition
[C{i},S{i}] = wavedec2(image{i},2,'bior3.7');
cfs2{i} = appcoef2(C{i},S{i},'bior3.7',2);
%gabor filter
fm{i} = featcal2(cfs2{i},5,8);
%gray-level co-occurence features extraction
offsets=[0 1;0 2;
    -1 1;-2 2;
    -1 0;-2 0;
    -1 -1;-2 -2];
glcms{i} = graycomatrix(cfs2{i},'offset',offsets);
stats{i} = graycoprops2(glcms{i},{'variance','entropy','IDM','Energy'});
stats2{i} = graycoprops(glcms{i},{'Contrast'});
cont{i}=stats2{i}.Contrast;
cont{i}=mean(cont{i});
v{i}=stats{i}.variance;
v{i}=mean(v{i});
IDM{i}=stats{i}.IDM;
IDM{i}=mean(IDM{i});
eng{i}=stats{i}.Energy;
eng{i}=mean(eng{i});
%gray-level run length features extraction
GLRLMS{i} = grayrlmatrix(cfs2{i},'NumLevels',255,'G',[0 256]);     
STATS{i} = grayrlprops(GLRLMS{i});       % rows=angles  colomns=properties
m{i}=mean(STATS{i});            %taking the mean of features  in each colomn
m{i}(:,1)=[];
m{i}(:,1)=[];
%features matrices
t(i,:)=cell2mat([eng(i) IDM(i) cont(i) fm(i) v(i) m(i)]);
end

