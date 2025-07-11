% used for the processing and plotting of MD simulation.
close all
clc
load("figure_4_colormaps.mat");

%% for angle 2, the horizontal position of the given data is good enough.
% try to plot it in the good way.
dat0=squeeze(s(71,:,:));        % 71 is the pore position.

dat0=permute(dat0,[2 1]);


% need to remove all the unnecessary data.
dat0(dat0<0.45)=nan;

max0=0.85;
dat0(dat0>max0)=max0;


% plot test.
imagesc(dat0);
axis xy equal;
colormap turbo


dat1(:,:,1)=dat0;
dat1(:,:,2)=dat0;


%% plotting part.


%% slicing for the water density.
%% rotate the matrix of density, angle 1. 51 degree difference between the two angles.
close all;
water_density_threshold_lower=0.003;
water_density_threshold_upper=0.07;

dx=imrotate3(density_theta2_nf,51,[0 0 1]);
slice=squeeze(dx(round(size(dx,1)/2),:,:));
slice(slice<water_density_threshold_lower)=nan;
slice(slice>water_density_threshold_upper)=water_density_threshold_upper;
h=imagesc(permute(slice,[2 1]));
set(h, 'AlphaData', ~isnan(permute(slice,[2 1])));
axis equal off xy

% colormap adjusted from 200 colormap MATLAB plugin, credits to Zhaoxu Liu/slandarer
colormap(water_density_blue);

%% rotate the matrix of density, angle 2. 51 degree difference between the two angles.
close all;
water_density_threshold_lower=0.05;
water_density_threshold_upper=0.3;

dx=imrotate3(density_theta2_nf,0,[0 0 1]);
slice=squeeze(dx(round(size(dx,1)/2),:,:));
slice(slice<water_density_threshold_lower)=nan;
slice(slice>water_density_threshold_upper)=water_density_threshold_upper;
h=imagesc(permute(slice,[2 1]));
set(h, 'AlphaData', ~isnan(permute(slice,[2 1])));
axis equal off xy

% colormap adjusted from 200 colormap MATLAB plugin, credits to Zhaoxu Liu/slandarer
colormap(water_density_blue);

%% full plot for the water orientation.
dx=imrotate3(orientation_theta2_nf,0,[0 0 1]);
slice=squeeze(dx(round(size(dx,1)/2),:,:));
% slice=smoothdata(slice,'gaussian',10);
h=imagesc(permute(slice,[2 1]));
set(h, 'AlphaData', ~isnan(permute(slice,[2 1])));
axis equal off xy
clim([-0.4 0.4]);

% colormap adjusted from 200 colormap MATLAB plugin, credits to Zhaoxu Liu/slandarer
colormap(water_orientation);

%% plotting the water dwell time, angle 1.
close all;
water_density_threshold_lower=0.01;
water_density_threshold_upper=0.3;

dx=imrotate3(data_theta2_nf,51,[0 0 1]);
dx=imgaussfilt3(dx,0.5);
slice=squeeze(dx(round(size(dx,1)/2),:,:));
slice(slice<water_density_threshold_lower)=nan;
slice(slice>water_density_threshold_upper)=water_density_threshold_upper;
h=imagesc(permute(slice,[2 1]));
set(h, 'AlphaData', ~isnan(permute(slice,[2 1])));
axis equal off xy

% colormap adjusted from 200 colormap MATLAB plugin, credits to Zhaoxu Liu/slandarer
colormap(water_density_blue);

%% plotting the water dwell time, angle 2.
close all;
water_density_threshold_lower=0.01;
water_density_threshold_upper=0.3;

dx=imrotate3(data_theta2_nf,0,[0 0 1]);
dx=imgaussfilt3(dx,0.5);
slice=squeeze(dx(round(size(dx,1)/2),:,:));
slice(slice<water_density_threshold_lower)=nan;
slice(slice>water_density_threshold_upper)=water_density_threshold_upper;
h=imagesc(permute(slice,[2 1]));
set(h, 'AlphaData', ~isnan(permute(slice,[2 1])));
axis equal off xy

% colormap adjusted from 200 colormap MATLAB plugin, credits to Zhaoxu Liu/slandarer
colormap(water_density_blue);


%% write tiff
tiff_write(data_theta2_nf,'D:\paper writing\3d afm\figure\MD simulation from Chris\2025-6-16 dwell time\1A_gaussian_filters_angles\1A.tiff')
%% locate the pore position, plot the density and the orientation. in 3D.
dx=orientation_theta2_nf;

% for a given threshold in density, it is viewed as no water density or
% part of the protein/lipid volume.
dx(density_theta2_nf<0.007)=nan;
dx=imrotate3(dx,45,[0 0 1]);
% cut a quarter of the data first for demonstration.
dx1=dx(1:round(size(dx,1)/2),1:round(size(dx,1)/2),:);
dx2=dx1(250:end,250:end,300:550);

%% another method, try to plot a cross section with the mean value.
% use the first dimension, all the first dimension falls into a single
% slice.
close all
dx1=dx2;


if 0        % use the max projection.
    dx=abs(dx1);
    [~,ind]=max(dx,[],1);


    [~, M, N] = size(dx1);
    [row, col] = ndgrid(1:M, 1:N);
    slice = dx1(sub2ind(size(dx1), squeeze(ind), row, col));
end

if 1
    % dx1(dx1==0)=nan;      % used to remove the likely no water region,
    % replaced with the comparison against water density matrix.
    slice=squeeze(mean(dx1,1,'omitmissing'));
end

slice=permute(slice,[2 1]);
% slice1=slice(300:550,250:end);

h=imagesc(slice);
colormap(water_orientation_narrow);
clim([-0.5 0.5]);
axis xy equal

set(h, 'AlphaData', ~isnan(slice));       % set the nan values to be transparent.


%% distribution of the orientation by layer

dx1=dx2;
close all

% try to do a auto sliced histogram on the z axis.

Nstart=1;
Nend=250;
interval=25;
hist=[];
N1=Nstart;

for i=1:ceil((Nend-Nstart)/interval)
    N2=N1+interval;
    if N2>Nend
        N2=Nend;
    end
    h=histogram(dx1(:,:,N1:N2));
    x=h.BinEdges;
    y=h.BinCounts;
    x=(x(1:end-1)+x(2:end))/2;
    hist{i}.x=x;
    hist{i}.y=y;
    N1=N2; 
end

figure;
t=tiledlayout(length(hist),1);
t.TileSpacing = 'none'; % 减少子图间距
t.Padding = 'compact';     % 减少整体边距

for i=length(hist):-1:1
    h(i)=nexttile;
    plot(hist{i}.x,hist{i}.y);
    % xline(0,'--',num2str(i));
end

linkaxes(h,'x');

for i=1:length(h)-1
    h(i).XAxis.Visible = 'on';
    h(i).Box = 'on';
end
xlabel('Water orientation');
ylabel('Counts (a.u.)');
