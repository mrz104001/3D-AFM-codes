% figure 2 script.


%% tiff reader.

clc;
clearvars -except path;

% must define the format of the file to be read.
file_format='.tif';

% read the file based on the previous loading and file format.
[name0,path0]=uigetfile({['*',file_format]},['Choose ',file_format,' file to read'],path);
if name0~=0
    path=[path0,name0];
else
    disp('File open error');
    return;
end

% double check if the file type is right or not.
if ~strcmp(name0(end-length(file_format)+1:end),file_format)
    disp(['File format is wrong, should be ',file_format,' files.']);
    return;
end

clear file_format;

% read the tiff file.

s0=tiffreadVolume(path);


%% do something to s0.
rotate_angle=-34;
s1=imrotate3(s0,rotate_angle,[0 0 1],"cubic","loose","FillValues",NaN);


%% give the coordinate for plotting a vertical slice in the data matrix.
j=47;
v_slice=permute(squeeze(s1(j,:,:)),[2 1]);
figure;
imagesc(v_slice);
colormap turbo
axis xy square off


%% added panel (hidden) for the plot of line profile across the AqpZ patches.
% the data is from C:\Runze\paper writing\3d afm\figures\Figure 2\AqpZ
% patches\height measurement\height across AqpZ patches.mat;

figure;
plot(height(:,1),height(:,2));
xlabel('x (nm)');
ylabel('height (nm)');
axis tight;