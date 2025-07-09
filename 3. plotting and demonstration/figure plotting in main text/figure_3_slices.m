% figure 3 making scripts.
clear
close all;

%% part I, input of the raw data, do an step average on z.

% the data will be used here is from test 5 calculation.
% file name would be
% structset_AqpZ_test5_0807a_v3_F2_force_20-Aug-2024-17-48-44-0400.mat
% the data used will be ppvolumeF2. the force data on a single AqpZ.

% do a step average on z direction, give a ppvolumeF3.

s=ppvolumeF3;

N_avg=5;            % the average step.
N_lowest=23;        % set the lowest point for full horizontal slice.
Ns=ceil(size(s,3)/N_avg);

s1=[];
for i=1:size(s,1)
    for j=1:size(s,2)
        for k=1:Ns
            s1(i,j,k)=mean(s(i,j,1+(k-1)*N_avg:min(k*N_avg,size(s,3))),'omitmissing');
        end
    end
end

ppvolumeF3=s1;

clearvars -except ppvolumeF3 N_avg;


%% get the force topography based on the largest slope of force.
% this will serve as the purpose of align the size and angle of data to
% protein structures.

s=ppvolumeF3;


for i=1:size(s,1)
    for j=1:size(s,2)
        [~,s_max0]=min(diff(s(i,j,:)));
        force_topo(i,j)=s_max0;
    end
end

clear s s_max0 i j;

% do a circle cut for the data to remove the outside feature.
cut_size=110;

close all
s_max2=force_topo;
s_c=size(s_max2,1)/2;
[xx,yy]=meshgrid(1:size(s_max2,1),1:size(s_max2,2));
s_max2(sqrt((xx-s_c).^2+(yy-s_c).^2)>cut_size)=nan;
force_topo=s_max2;

clear s_c xx yy s_max2

slice_h=force_topo;

% try to remove the nan values on the side.
slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
slice_h(:,isnan(slice_h_temp))=[];

slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
slice_h(isnan(slice_h_temp),:)=[];

imagesc(slice_h);
axis xy equal off
colormap turbo

clear slice_h slice_h_temp

% optional, copy the graph to the system.
% copygraphics(gcf, 'BackgroundColor', 'none', 'ContentType', 'image','Resolution',300);


%% get an isosurface of force plot at the bottom.

% get the largest value for each pixel.

s = ppvolumeF3;

for i=1:size(s,1)
    for j=1:size(s,2)
        s_max(i,j)=max(s(i,j,:));
    end
end

s_max(s_max==0)=nan;
% s_max(s_max<0.3*mean(s_max,'all','omitmissing'))=nan;
s_max0=0.3*mean(s_max,'all','omitmissing');


for i=1:size(s,1)
    for j=1:size(s,2)
        isoforce(i,j)=find(s(i,j,:)<s_max0,1,'first');
    end
end
isoforce(isnan(s_max))=nan;

cut_size=110;

close all
s_max2=isoforce;
s_c=size(s_max2,1)/2;
[xx,yy]=meshgrid(1:size(s_max2,1),1:size(s_max2,2));
s_max2(sqrt((xx-s_c).^2+(yy-s_c).^2)>cut_size)=nan;
isoforce=s_max2;


% remove nan values on the sides.
s_max1_temp=isoforce(ceil(size(isoforce,1)/2),:);
isoforce(:,isnan(s_max1_temp))=[];

s_max1_temp=isoforce(:,ceil(size(isoforce,2)/2));
isoforce(isnan(s_max1_temp),:)=[];

imagesc(isoforce);
axis xy equal off
colormap turbo

clear i j s_max s s_c s_max0 s_max1_temp s_max2 xx yy

%% get the given height of the data.
close all
N_horizontal=18;

s=ppvolumeF3;
slice_h=squeeze(s(:,:,N_horizontal));

% try to remove the nan values on the side.
slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
slice_h(:,isnan(slice_h_temp))=[];

slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
slice_h(isnan(slice_h_temp),:)=[];

imagesc(slice_h);
axis xy equal off
colormap turbo

clear slice_h

%% do an image rotation. only for vertical slices. rotate 2: cut through pore channel.
close all
rotate_angle_2=27.26;            % changed from 28 to 22, then to 26, more on the pore.
ppvolumeF3_v=imrotate3(ppvolumeF3,rotate_angle_2,[0 0 1],"cubic","loose","FillValues",NaN);
slice_h=squeeze(ppvolumeF3_v(:,:,N_lowest));

% try to remove the nan values on the side.
slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
slice_h(:,isnan(slice_h_temp))=[];

slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
slice_h(isnan(slice_h_temp),:)=[];

imagesc(slice_h);
axis xy equal off
colormap turbo

cmin=min(ppvolumeF3_v,[],'all');
cmax=max(ppvolumeF3_v,[],'all');

figure;
slice_v=permute(squeeze(ppvolumeF3_v(ceil(size(ppvolumeF3_v,1)/2),:,:)),[2 1]);

% try to remove all the nan on the side.
slice_v_temp=slice_v(end,:);
slice_v(:,isnan(slice_v_temp))=[];
h=imagesc(slice_v);
axis xy off
daspect([8 1 1])           % roughly set the ratio to 2:1 for the plot in Illustrator.
clim manual
clim([cmin cmax]);
colormap turbo

set(h, 'AlphaData', ~isnan(slice_v));       % set the nan values to be transparent.


% clear all the unwanted variables.
clear h slice_h slice_h_temp slice_v slice_v_temp ppvolumeF3_v

%% do an image rotation. only for vertical slices. rotate 1: cut through protein peak.
close all
rotate_angle_1=-19.74;           % changed from -20 to -29 for better alignment on highest point.
ppvolumeF3_v=imrotate3(ppvolumeF3,rotate_angle_1,[0 0 1],"cubic","loose","FillValues",NaN);
slice_h=squeeze(ppvolumeF3_v(:,:,N_lowest));

% try to remove the nan values on the side.
slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
slice_h(:,isnan(slice_h_temp))=[];

slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
slice_h(isnan(slice_h_temp),:)=[];

imagesc(slice_h);
axis xy equal off
colormap turbo


cmin=min(ppvolumeF3_v,[],'all');
cmax=max(ppvolumeF3_v,[],'all');

figure;
slice_v=permute(squeeze(ppvolumeF3_v(ceil(size(ppvolumeF3_v,1)/2),:,:)),[2 1]);

% try to remove all the nan on the side.
slice_v_temp=slice_v(end,:);
slice_v(:,isnan(slice_v_temp))=[];
h=imagesc(slice_v);
axis xy off
daspect([8 1 1])           % roughly set the ratio to 2:1 for the plot in Illustrator.
clim manual
clim([cmin cmax]);
colormap turbo

set(h, 'AlphaData', ~isnan(slice_v));       % set the nan values to be transparent.

% clear all the unwanted variables.
clear h slice_h slice_h_temp slice_v slice_v_temp ppvolumeF3_v

%% do an image rotation. only for vertical slices. rotate 3: cut through protein peak.
close all
rotate_angle_3=0;
ppvolumeF3_v=imrotate3(ppvolumeF3,rotate_angle_3,[0 0 1],"cubic","loose","FillValues",NaN);
slice_h=squeeze(ppvolumeF3_v(:,:,N_lowest));

% try to remove the nan values on the side.
slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
slice_h(:,isnan(slice_h_temp))=[];

slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
slice_h(isnan(slice_h_temp),:)=[];

imagesc(slice_h);
axis xy equal off
colormap turbo


cmin=min(ppvolumeF3_v,[],'all');
cmax=max(ppvolumeF3_v,[],'all');

figure;
slice_v=permute(squeeze(ppvolumeF3_v(ceil(size(ppvolumeF3_v,1)/2),:,:)),[2 1]);

% try to remove all the nan on the side.
slice_v_temp=slice_v(end,:);
slice_v(:,isnan(slice_v_temp))=[];
h=imagesc(slice_v);
axis xy off
daspect([8 1 1])           % roughly set the ratio to 2:1 for the plot in Illustrator.
clim manual
clim([cmin cmax]);
colormap turbo

set(h, 'AlphaData', ~isnan(slice_v));       % set the nan values to be transparent.

% clear all the unwanted variables.
clear h slice_h slice_h_temp slice_v slice_v_temp ppvolumeF3_v

%% use the original force matrix and save all the height slices given in a coordinate.
close all;
% h=[23,28,37,45,53,41,50];

% h=[23 31 35 40 45];         % new height slices.

% h=[23 31 36 42 48];

h=[21 31 33 36 42 48];         % 2025.5.2, use the lowest point as z = 0.

% for the given original data, give the same color limits.
cmin=min(ppvolumeF3(:,:,h(1)),[],'all');
cmax=max(ppvolumeF3(:,:,h(1)),[],'all');

for i=1:length(h)
    figure;
    slice_h=squeeze(ppvolumeF3(:,:,h(i)));
    % remove nan values on the sides.
    slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
    slice_h(:,isnan(slice_h_temp))=[];

    slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
    slice_h(isnan(slice_h_temp),:)=[];

    imagesc(slice_h);
    colormap turbo
    clim manual
    clim([cmin cmax]);
    axis xy equal off
end

clear i slice_h slice_h_temp

%% help with calculation of the actual height in A reference to the lowest z.
z_res=zres_tgt*N_avg;
hz=(h-h(1)).*z_res;

% %% get the line profile for the simulated PDB surface of AqpZ. not in use currently.
% 
% % the files are saved in C:\Runze\paper writing\3d afm\figures\Figure 3\3d
% % LAFM\PDB generated surface\measures.
% % load the file before plotting it.
% 
% x=lineprofile(:,1);
% h=lineprofile(:,2);
% 
% figure;
% plot(x/10,h/10);
% 
% axis tight equal
% xlabel('X (nm)');
% ylabel('z (nm)');
% hold on;
% 
% % try to do a tip convolution and see what that is.
% % the basic idea is to use a circle to find the tangent point.
% 
% 
% tip_r=50;            % tip radius in angstrom.
% N=100;
% 
% h_conv=zeros(length(h),1);
% h_max=max(h);       % highest point in the data. in angstrom.
% h_min=min(h);       % lowest point in the data.
% h_ini=h_max+tip_r+10;    % all the search starts from this point.
% h_end=h_min+tip_r;      % all the search ends before this point.
% 
% for i=1:length(h)
%     for j=1:N
%         h_t=h_ini+(h_end-h_ini)/(N-1)*(j-1);
% 
%         % calculate the distance of every point in the curve to the tip.
%         h_d=sqrt((x-i).^2+(h-h_t).^2);
%         if min(h_d)<=tip_r
%             h_conv(i)=h_t-tip_r;
%             break;
%         end
%     end
% end
% 
% plot(x/10,h_conv/10,'-');
% hold off;

%% change colormap for the good plot.   vertical slices.
colormap(water_density_blue);
clim([-1e-11 1e-11]);


%% change colormap for the good plot.   horizontal slices.
colormap(water_density_blue);
clim([-2e-12 1e-11]);



%% horizontal slices. auto
colormap(water_density_blue);
clim auto;


%% for panel c, add the local minimum position.
close all;

% for the given height vector, do plotting and find the local minimum.
% the height starts from 21. every 3 pixels.
h=[21 26 31 36 42 48];         % 2025.5.2, use the lowest point as z = 0.

% h=[36];

% for the given original data, give the same color limits.
cmin=min(ppvolumeF3(:,:,h(1)),[],'all');
cmax=max(ppvolumeF3(:,:,h(1)),[],'all');

local_min={};

for i=1:length(h)
    figure('Name',['h=',num2str(h(i))]);
    slice_h=squeeze(ppvolumeF3(:,:,h(i)));
    % remove nan values on the sides.
    slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
    slice_h(:,isnan(slice_h_temp))=[];

    slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
    slice_h(isnan(slice_h_temp),:)=[];

    
    tf=islocalmin2(slice_h,"MinSeparation",10,'MinProminence',0.1e-12);

    % to remove the features that are far away from the center. 85 is
    % determined by the size of the AqpZ, only inside the AqpZ excellular
    % side is included in this data.
    tf=processMatrix(tf,70,0);

    % find all the coordinates of the local minima.
    [Y,X]=find(tf);

    % find the minimum of it.
    z=[];
    for zi=1:length(Y)
        z(zi)=slice_h(Y(zi),X(zi));
    end

    % find the lowest local minima as the global minima. with tolerance.
    ind=find(ismembertol(z,min(z),1e-12));

    local_min{end+1}.X=X;
    local_min{end}.Y=Y;
    local_min{end}.h=h(i);
    
    imagesc(slice_h);
    colormap(water_density_blue);
    clim auto;
    axis xy equal off

    hold on;
    for pi=1:length(X)
        if any(ismember(ind,pi),'all')
            plot(X(pi),Y(pi),'x','Color','r');
        else
            plot(X(pi),Y(pi),'x','Color','b');
        end
    end
    hold off;
end

% clear i slice_h slice_h_temp

%% for panel c, corresponding SI figure.
close all;

% for the given height vector, do plotting and find the local minimum.
% the height starts from 21. every 3 pixels.
h=21:48;         % 2025.5.2, use the lowest point as z = 0.

z_res=zres_tgt*N_avg;
hz=(h-h(1)).*z_res;

% create a tilelayout.
fig=figure;
tiledlayout(4,7,'TileSpacing','tight','Padding','tight');



local_min={};
cm={};

for i=1:length(h)
    nexttile;
    slice_h=squeeze(ppvolumeF3(:,:,h(i)));
    % remove nan values on the sides.
    slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
    slice_h(:,isnan(slice_h_temp))=[];

    slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
    slice_h(isnan(slice_h_temp),:)=[];

    
    tf=islocalmin2(slice_h,"MinSeparation",10,'MinProminence',0.01e-12);

    % to remove the features that are far away from the center. 85 is
    % determined by the size of the AqpZ, only inside the AqpZ excellular
    % side is included in this data.
    tf=processMatrix(tf,70,0);

    % find all the coordinates of the local minima.
    [Y,X]=find(tf);

    % find the minimum of it.
    z=[];
    for zi=1:length(Y)
        z(zi)=slice_h(Y(zi),X(zi));
    end
    z=z-max(z);

    % find the lowest local minima as the global minima. with tolerance.
    
    % z_threshold=0.5;
    % ind1=find(z>min(z)*z_threshold);
    % X(ind1)=[];
    % Y(ind1)=[];
    % z(ind1)=[];

    ind=find(ismembertol(z,min(z),1e-12));
    local_min{end+1}.X=X;
    local_min{end}.Y=Y;
    local_min{end}.h=h(i);

    cm{end+1}.min=min(slice_h,[],'all');
    cm{end}.max=max(slice_h,[],'all');
    
    imagesc(slice_h);
    colormap(water_density_blue);
    clim auto;
    axis xy equal off

    hold on;
    for pi=1:length(X)
        if any(ismember(ind,pi),'all')
            plot(X(pi),Y(pi),'x','Color','r');
        else
            plot(X(pi),Y(pi),'x','Color','b');
        end
    end
    hold off;
    title(['z = ',num2str(hz(i),2)]);
end

exportgraphics(fig, 'full slices at all heights.pdf', 'ContentType', 'vector','Resolution',600);

% clear i slice_h slice_h_temp


%% choose the right heights for the panel b.
z=[21 27 34 36 37 41 44 48];

% for the side plot, try to line up.
close all
rotate_angle_2=27.26;            % changed from 28 to 22, then to 26, more on the pore.
ppvolumeF3_v=imrotate3(ppvolumeF3,rotate_angle_2,[0 0 1],"cubic","loose","FillValues",NaN);

figure;
slice_v=permute(squeeze(ppvolumeF3_v(ceil(size(ppvolumeF3_v,1)/2),:,:)),[2 1]);

% try to remove all the nan on the side.
slice_v_temp=slice_v(end,:);
slice_v(:,isnan(slice_v_temp))=[];
h=imagesc(slice_v);
axis xy off
daspect([8 1 1])           % roughly set the ratio to 2:1 for the plot in Illustrator.
clim manual
clim([cmin cmax]);
colormap turbo

hold on
yline(z,'--w');
hold off

set(h, 'AlphaData', ~isnan(slice_v));       % set the nan values to be transparent.


% clear all the unwanted variables.
clear h slice_h slice_h_temp slice_v slice_v_temp ppvolumeF3_v

%% for panel d, have all the local minima and record it in a cell variable.
% also gives a movies as SI.

% still working on gettting the right dots as local or global minima.

close all;

% for the given height vector, do plotting and find the local minimum.
% the height starts from 21. every 3 pixels.
h=21:1:48;

z_res=zres_tgt*N_avg;
hz=(h-h(1)).*z_res;
hz=roundn(hz,-1);

% h=[36];

% for the given original data, give the same color limits.
cmin=min(ppvolumeF3(:,:,h(1)),[],'all');
cmax=max(ppvolumeF3(:,:,h(1)),[],'all');

local_min={};

% define the movie settings.
v = VideoWriter('local_minima.mp4', 'MPEG-4');
v.FrameRate = 10;  % Adjust as needed
open(v);

% Create output folder for TIFFs
tiff_filename = 'local_minima.tiff';

figure_handle=figure('Name',['h=',num2str(h(i))]);

for i=1:length(h)
    slice_h=squeeze(ppvolumeF3(:,:,h(i)));
    % remove nan values on the sides.
    slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
    slice_h(:,isnan(slice_h_temp))=[];

    slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
    slice_h(isnan(slice_h_temp),:)=[];

    
    tf=islocalmin2(slice_h,"MinSeparation",5,'MinProminence',0.01e-12);

    % to remove the features that are far away from the center. 85 is
    % determined by the size of the AqpZ, only inside the AqpZ excellular
    % side is included in this data.
    tf=processMatrix(tf,70,0);

    % find all the coordinates of the local minima.
    [Y,X]=find(tf);

    % find the minimum of it.
    z=[];
    for zi=1:length(Y)
        z(zi)=slice_h(Y(zi),X(zi));
    end

    % find the lowest local minima as the global minima. with tolerance.
    ind=find(ismembertol(z,min(z),1e-12));

    local_min{end+1}.X=X;
    local_min{end}.Y=Y;
    local_min{end}.h=h(i);
    local_min{end}.ind=ind;
    
    imagesc(slice_h);
    colormap(water_density_blue);
    clim auto;
    axis xy equal off

    hold on;
    for pi=1:length(X)
        if any(ismember(ind,pi),'all')
            h1=plot(X(pi),Y(pi),'x','Color','r');

            % change the marker size and thickness of global minima.
            h1.LineWidth=1.5;
            h1.MarkerSize=7;
        else
            h2=plot(X(pi),Y(pi),'x','Color','b');
        end
    end
    hold off;

    % set the title for each frame.
    title(['h = ',num2str(hz(i),'%0.1f'),' Å']);

    drawnow;           % Update the figure
    % Record current frame
    frame = getframe(figure_handle);
    writeVideo(v, frame);

    % Convert frame to image
    [im, ~] = frame2im(frame);
    
    % Save as multi-page TIFF
    if i == 1
        imwrite(im, tiff_filename, 'tiff', 'Compression', 'none');
    else
        imwrite(im, tiff_filename, 'tiff', 'WriteMode', 'append', 'Compression', 'none');
    end
end

close(v);

% clear i slice_h slice_h_temp


%% plot the above into a 3D demo of all the minima.

% first convert the cell variable to a [x,y,z] matrix.
close all;
pos=[];
for i=1:length(local_min)
    pos=[pos;local_min{i}.X,local_min{i}.Y,local_min{i}.h.*ones(length(local_min{i}.X),1)];
end

% connect all the global minima index into one vector.
pos_GlobalMin=[];
pos_temp=0;
for i=1:length(local_min)
    if i>1
        pos_temp=pos_temp+length(local_min{i-1}.X);
    end
    pos_GlobalMin=[pos_GlobalMin,local_min{i}.ind+pos_temp];
end

% 
groups=connectLayerPaths5(pos,pos_GlobalMin,20);
groups_sub=divideGroupsFourFold(pos,groups);

daspect([8 8 1]);

% write it to the local folder.
fig = gcf;  % or specify your figure handle
exportgraphics(fig, 'Connected_local_minima.pdf', 'ContentType', 'vector');






%% test of the subgroups, correct.
close all

figure;
hold on

j=2;

for i=1:length(groups_sub{j})
    ind=groups_sub{j}{i};
    if length(ind)>2
        plot3(pos(ind,1),pos(ind,2),pos(ind,3),'o-');
    end
end


daspect([8 8 1]);



%% fit the local minima.
close all;
peak=[];
R=[];
for i=1:size(pos,1)
    slice_h=squeeze(density(:,:,pos(i,3)-20));
    % remove nan values on the sides.
    slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
    slice_h(:,isnan(slice_h_temp))=[];

    slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
    slice_h(isnan(slice_h_temp),:)=[];

    widsize = 3;
    rad = widsize + 1;
    while rad > widsize - 0.1
        widsize = widsize + 2;
        [peak(i),R(i)]=fitSymmetricGaussian2DAtPoint(slice_h,[pos(i,1),pos(i,2)],widsize, 1);
        rad = peak(i);
    end
end


%% check on the local minima.
close all
% Create output folder for TIFFs
tiff_filename = 'peak_fit.tiff';
figure_handle=figure;
for i=1:length(peak)
    slice_h=squeeze(density(:,:,pos(i,3)-20));
    slice_h_temp=slice_h(ceil(size(slice_h,1)/2),:);
    slice_h(:,isnan(slice_h_temp))=[];

    slice_h_temp=slice_h(:,ceil(size(slice_h,2)/2));
    slice_h(isnan(slice_h_temp),:)=[];

    imagesc(slice_h);
    axis xy off equal
    hold on
    plot(pos(i,1),pos(i,2),'x','Color','r');

    center = [pos(i,1) pos(i,2)];
    radius = peak(i);
    x = center(1) - radius;
    y = center(2) - radius;
    w = radius * 2;
    h = radius * 2;
    pos1 = [x y w h];
    cur = [1 1];

    rectangle('Position',pos1,'Curvature',cur)

    hold off

    frame = getframe(figure_handle);
    % Convert frame to image
    [im, ~] = frame2im(frame);
    
    % Save as multi-page TIFF
    if i == 1
        imwrite(im, tiff_filename, 'tiff', 'Compression', 'none');
    else
        imwrite(im, tiff_filename, 'tiff', 'WriteMode', 'append', 'Compression', 'none');
    end
end

%% local functions.
function A_processed = processMatrix(A, R, s)
    % Get the size of the matrix
    [rows, cols] = size(A);
    
    % Compute the center coordinates (can be fractional for even sizes)
    center_row = (rows + 1) / 2;
    center_col = (cols + 1) / 2;
    
    % Create a meshgrid of coordinates
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Compute the distance of each point from the center
    dist_from_center = sqrt((X - center_col).^2 + (Y - center_row).^2);
    
    % Create a mask for elements farther than R from the center
    mask = dist_from_center > R;
    
    % Copy the original matrix
    A_processed = A;
    
    % Apply the value s to elements outside the radius R
    A_processed(mask) = s;
end

function connectLayerPaths(P)
    % Sort by Z (height)
    % P = sortrows(P, 3); % since it is already in order, neglect this
    % step.
    layers = unique(P(:,3));
    n = size(P,1);

    % Build point index
    pointIdx = (1:n)';
    P = [P, pointIdx];  % Append index as column 4

    % Initialize connection map: child_idx -> parent_idx(s)
    connections = containers.Map('KeyType','int32','ValueType','any');

    for k = 2:length(layers)
        z_current = layers(k);
        z_lower = layers(k-1);
        current_layer = P(P(:,3) == z_current, :);
        lower_layer = P(P(:,3) == z_lower, :);

        for i = 1:size(current_layer,1)
            pt = current_layer(i,1:2);
            distances = vecnorm(lower_layer(:,1:2) - pt, 2, 2);
            min_dist = min(distances);
            tol = 2;
            closest_idxs = find(abs(distances - min_dist) < tol);

            for j = 1:length(closest_idxs)
                child_id = current_layer(i,4);
                parent_id = lower_layer(closest_idxs(j),4);
                if isKey(connections, child_id)
                    connections(child_id) = [connections(child_id), parent_id];
                else
                    connections(child_id) = parent_id;
                end
            end
        end
    end

    % Find root nodes (points in bottom layer)
    root_ids = P(P(:,3) == layers(1), 4);
    cmap = lines(length(root_ids));
    paths = containers.Map('KeyType','int32','ValueType','any');

    % Assign a unique color to each root path
    for i = 1:length(root_ids)
        queue = root_ids(i);
        paths(root_ids(i)) = cmap(i,:);
        while ~isempty(queue)
            curr = queue(1);
            queue(1) = [];
            for child = keys(connections)
                cid = child{1};
                parents = connections(cid);
                if any(parents == curr) && ~isKey(paths, cid)
                    paths(cid) = paths(curr);
                    queue(end+1) = cid;
                end
            end
        end
    end

    % Plot
    figure;
    hold on;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    scatter3(P(:,1), P(:,2), P(:,3), 36, 'k', 'filled');

    % Plot lines with path colors
    for cid = keys(connections)
        child_id = cid{1};
        parent_ids = connections(child_id);
        child_point = P(P(:,4)==child_id, 1:3);
        for j = 1:length(parent_ids)
            parent_point = P(P(:,4)==parent_ids(j), 1:3);
            clr = paths(child_id);
            plot3([child_point(1), parent_point(1)], ...
                  [child_point(2), parent_point(2)], ...
                  [child_point(3), parent_point(3)], ...
                  '-', 'Color', clr, 'LineWidth', 2);
        end
    end

    hold off;
end


function connectLayerPaths2(P, pos_GlobalMin)
    % Sort by Z (height)
    layers = unique(P(:,3));
    n = size(P,1);

    % Build point index
    pointIdx = (1:n)';
    P = [P, pointIdx];  % Append index as column 4

    % Initialize connection map: child_idx -> parent_idx(s)
    connections = containers.Map('KeyType','int32','ValueType','any');

    for k = 2:length(layers)
        z_current = layers(k);
        z_lower = layers(k-1);
        current_layer = P(P(:,3) == z_current, :);
        lower_layer = P(P(:,3) == z_lower, :);

        for i = 1:size(current_layer,1)
            pt = current_layer(i,1:2);
            distances = vecnorm(lower_layer(:,1:2) - pt, 2, 2);
            min_dist = min(distances);
            tol = 2;
            closest_idxs = find(abs(distances - min_dist) < tol);

            for j = 1:length(closest_idxs)
                child_id = current_layer(i,4);
                parent_id = lower_layer(closest_idxs(j),4);
                if isKey(connections, child_id)
                    connections(child_id) = [connections(child_id), parent_id];
                else
                    connections(child_id) = parent_id;
                end
            end
        end
    end

    % Find root nodes (points in bottom layer)
    root_ids = P(P(:,3) == layers(1), 4);
    cmap = lines(length(root_ids));
    paths = containers.Map('KeyType','int32','ValueType','any');

    % Assign a unique color to each root path
    for i = 1:length(root_ids)
        queue = root_ids(i);
        paths(root_ids(i)) = cmap(i,:);
        while ~isempty(queue)
            curr = queue(1);
            queue(1) = [];
            for child = keys(connections)
                cid = child{1};
                parents = connections(cid);
                if any(parents == curr) && ~isKey(paths, cid)
                    paths(cid) = paths(curr);
                    queue(end+1) = cid;
                end
            end
        end
    end

    % Plot
    figure;
    hold on;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % Highlighted points (larger size)
    scatter3(P(pos_GlobalMin,1), P(pos_GlobalMin,2), P(pos_GlobalMin,3), ...
             100, 'r', 'filled');  % Large red markers for GlobalMin

    % All other points (normal size)
    other_idx = setdiff(1:n, pos_GlobalMin);
    scatter3(P(other_idx,1), P(other_idx,2), P(other_idx,3), ...
             36, 'k', 'filled');  % Small black markers

    % Plot lines with path colors
    for cid = keys(connections)
        child_id = cid{1};
        parent_ids = connections(child_id);
        child_point = P(P(:,4)==child_id, 1:3);
        for j = 1:length(parent_ids)
            parent_point = P(P(:,4)==parent_ids(j), 1:3);
            clr = paths(child_id);
            plot3([child_point(1), parent_point(1)], ...
                  [child_point(2), parent_point(2)], ...
                  [child_point(3), parent_point(3)], ...
                  '-', 'Color', clr, 'LineWidth', 2);
        end
    end

    hold off;
end


function connectLayerPaths4(P, pos_GlobalMin, dist_thresh)
    % Sort by Z (height)
    layers = unique(P(:,3));
    n = size(P,1);

    % Build point index
    pointIdx = (1:n)';
    P = [P, pointIdx];  % Append index as column 4

    % Initialize connection map: child_idx -> parent_idx(s)
    connections = containers.Map('KeyType','int32','ValueType','any');

    for k = 2:length(layers)
        z_current = layers(k);
        z_lower = layers(k-1);
        current_layer = P(P(:,3) == z_current, :);
        lower_layer = P(P(:,3) == z_lower, :);

        for i = 1:size(current_layer,1)
            pt = current_layer(i,1:2);
            distances = vecnorm(lower_layer(:,1:2) - pt, 2, 2);
            min_dist = min(distances);

            % Skip if the closest point is beyond the threshold
            if min_dist > dist_thresh
                continue;
            end

            tol = 2;
            closest_idxs = find(abs(distances - min_dist) < tol);

            for j = 1:length(closest_idxs)
                child_id = current_layer(i,4);
                parent_id = lower_layer(closest_idxs(j),4);
                if isKey(connections, child_id)
                    connections(child_id) = [connections(child_id), parent_id];
                else
                    connections(child_id) = parent_id;
                end
            end
        end
    end

    % Assign a unique color to every connected component
    visited = false(n,1);
    paths = containers.Map('KeyType','int32','ValueType','any');
    cmap = lines(n);  % Enough colors for all points

    color_idx = 1;
    for i = 1:n
        if visited(i)
            continue;
        end

        queue = i;
        paths(i) = cmap(color_idx,:);
        visited(i) = true;

        while ~isempty(queue)
            curr = queue(1);
            queue(1) = [];

            % Check children: who have curr as parent
            for child = keys(connections)
                cid = child{1};
                parents = connections(cid);
                if any(parents == curr) && ~isKey(paths, cid)
                    paths(cid) = cmap(color_idx,:);
                    queue(end+1) = cid;
                    visited(cid) = true;
                end
            end

            % Check parents of curr too (for completeness)
            if isKey(connections, curr)
                parents = connections(curr);
                for p = parents
                    if ~visited(p)
                        paths(p) = cmap(color_idx,:);
                        queue(end+1) = p;
                        visited(p) = true;
                    end
                end
            end
        end

        color_idx = color_idx + 1;
    end

    % Plot
    figure;
    hold on;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % Highlighted points (larger size)
    scatter3(P(pos_GlobalMin,1), P(pos_GlobalMin,2), P(pos_GlobalMin,3), ...
             50, 'r', 'filled');  % Large red markers for GlobalMin

    % All other points (normal size)
    other_idx = setdiff(1:n, pos_GlobalMin);
    scatter3(P(other_idx,1), P(other_idx,2), P(other_idx,3), ...
             20, 'k', 'filled');  % Small black markers

    % Plot lines with path colors
    for cid = keys(connections)
        child_id = cid{1};
        parent_ids = connections(child_id);
        child_point = P(P(:,4)==child_id, 1:3);
        for j = 1:length(parent_ids)
            parent_point = P(P(:,4)==parent_ids(j), 1:3);
            if isKey(paths, child_id)
                clr = paths(child_id);
                plot3([child_point(1), parent_point(1)], ...
                      [child_point(2), parent_point(2)], ...
                      [child_point(3), parent_point(3)], ...
                      '-', 'Color', clr, 'LineWidth', 2);
            end
        end
    end

    hold off;
end


function groups = connectLayerPaths5(P, pos_GlobalMin, dist_thresh)
    layers = unique(P(:,3));
    n = size(P,1);
    pointIdx = (1:n)';
    P = [P, pointIdx];  % Append index as column 4

    connections = containers.Map('KeyType','int32','ValueType','any');

    for k = 2:length(layers)
        z_current = layers(k);
        z_lower = layers(k-1);
        current_layer = P(P(:,3) == z_current, :);
        lower_layer = P(P(:,3) == z_lower, :);

        for i = 1:size(current_layer,1)
            pt = current_layer(i,1:2);
            distances = vecnorm(lower_layer(:,1:2) - pt, 2, 2);
            min_dist = min(distances);

            if min_dist > dist_thresh
                continue;
            end

            tol = 2;
            closest_idxs = find(abs(distances - min_dist) < tol);

            for j = 1:length(closest_idxs)
                child_id = current_layer(i,4);
                parent_id = lower_layer(closest_idxs(j),4);
                if isKey(connections, child_id)
                    connections(child_id) = [connections(child_id), parent_id];
                else
                    connections(child_id) = parent_id;
                end
            end
        end
    end

    % Find connected groups and assign colors
    visited = false(n,1);
    paths = containers.Map('KeyType','int32','ValueType','any');
    raw_groups = {};
    group_lengths = [];
    color_idx = 1;

    for i = 1:n
        if visited(i)
            continue;
        end

        queue = i;
        group = i;
        paths(i) = color_idx;
        visited(i) = true;

        while ~isempty(queue)
            curr = queue(1);
            queue(1) = [];

            for child = keys(connections)
                cid = child{1};
                parents = connections(cid);
                if any(parents == curr) && ~visited(cid)
                    paths(cid) = color_idx;
                    queue(end+1) = cid;
                    visited(cid) = true;
                    group(end+1) = cid;
                end
            end

            if isKey(connections, curr)
                parents = connections(curr);
                for p = parents
                    if ~visited(p)
                        paths(p) = color_idx;
                        queue(end+1) = p;
                        visited(p) = true;
                        group(end+1) = p;
                    end
                end
            end
        end

        % Compute path length of this group
        length_sum = 0;
        for j = 1:length(group)
            pt_id = group(j);
            if isKey(connections, pt_id)
                parents = connections(pt_id);
                pt_pos = P(P(:,4)==pt_id,1:3);
                for p = parents
                    p_pos = P(P(:,4)==p,1:3);
                    length_sum = length_sum + norm(pt_pos - p_pos);
                end
            end
        end

        raw_groups{color_idx} = group; %#ok<AGROW>
        group_lengths(color_idx) = length_sum; %#ok<AGROW>
        color_idx = color_idx + 1;
    end

    % Sort groups by descending path length
    [~, sorted_idx] = sort(group_lengths, 'descend');
    groups = raw_groups(sorted_idx);

    % Remap paths with new sorted color indices
    cmap = lines(length(groups));
    path_color_map = containers.Map('KeyType','int32','ValueType','any');
    for i = 1:length(groups)
        for pt_id = groups{i}
            path_color_map(pt_id) = cmap(i,:);
        end
    end

    % Plot
    figure;
    hold on;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    scatter3(P(pos_GlobalMin,1), P(pos_GlobalMin,2), P(pos_GlobalMin,3), ...
             100, 'r', 'filled');  % Highlighted markers

    other_idx = setdiff(1:n, pos_GlobalMin);
    scatter3(P(other_idx,1), P(other_idx,2), P(other_idx,3), ...
             36, 'k', 'filled');
    

    for cid = keys(connections)
        child_id = cid{1};
        parent_ids = connections(child_id);
        child_point = P(P(:,4)==child_id, 1:3);
        for j = 1:length(parent_ids)
            parent_point = P(P(:,4)==parent_ids(j), 1:3);
            if isKey(path_color_map, child_id)
                clr = path_color_map(child_id);
                plot3([child_point(1), parent_point(1)], ...
                      [child_point(2), parent_point(2)], ...
                      [child_point(3), parent_point(3)], ...
                      '-', 'Color', clr, 'LineWidth', 2);
            end
        end
    end

    hold off;
end

function groups = connectLayerPaths6(P, pos_GlobalMin, dist_thresh, peak)
    layers = unique(P(:,3));
    n = size(P,1);
    pointIdx = (1:n)';
    P = [P, pointIdx];  % Append index as column 4

    connections = containers.Map('KeyType','int32','ValueType','any');

    for k = 2:length(layers)
        z_current = layers(k);
        z_lower = layers(k-1);
        current_layer = P(P(:,3) == z_current, :);
        lower_layer = P(P(:,3) == z_lower, :);

        for i = 1:size(current_layer,1)
            pt = current_layer(i,1:2);
            distances = vecnorm(lower_layer(:,1:2) - pt, 2, 2);
            min_dist = min(distances);

            if min_dist > dist_thresh
                continue;
            end

            tol = 2;
            closest_idxs = find(abs(distances - min_dist) < tol);

            for j = 1:length(closest_idxs)
                child_id = current_layer(i,4);
                parent_id = lower_layer(closest_idxs(j),4);
                if isKey(connections, child_id)
                    connections(child_id) = [connections(child_id), parent_id];
                else
                    connections(child_id) = parent_id;
                end
            end
        end
    end

    % Find connected groups and assign colors
    visited = false(n,1);
    paths = containers.Map('KeyType','int32','ValueType','any');
    raw_groups = {};
    group_lengths = [];
    color_idx = 1;

    for i = 1:n
        if visited(i)
            continue;
        end

        queue = i;
        group = i;
        paths(i) = color_idx;
        visited(i) = true;

        while ~isempty(queue)
            curr = queue(1);
            queue(1) = [];

            for child = keys(connections)
                cid = child{1};
                parents = connections(cid);
                if any(parents == curr) && ~visited(cid)
                    paths(cid) = color_idx;
                    queue(end+1) = cid;
                    visited(cid) = true;
                    group(end+1) = cid;
                end
            end

            if isKey(connections, curr)
                parents = connections(curr);
                for p = parents
                    if ~visited(p)
                        paths(p) = color_idx;
                        queue(end+1) = p;
                        visited(p) = true;
                        group(end+1) = p;
                    end
                end
            end
        end

        % Compute path length of this group
        length_sum = 0;
        for j = 1:length(group)
            pt_id = group(j);
            if isKey(connections, pt_id)
                parents = connections(pt_id);
                pt_pos = P(P(:,4)==pt_id,1:3);
                for p = parents
                    p_pos = P(P(:,4)==p,1:3);
                    length_sum = length_sum + norm(pt_pos - p_pos);
                end
            end
        end

        raw_groups{color_idx} = group; %#ok<AGROW>
        group_lengths(color_idx) = length_sum; %#ok<AGROW>
        color_idx = color_idx + 1;
    end

    % Sort groups by descending path length
    [~, sorted_idx] = sort(group_lengths, 'descend');
    groups = raw_groups(sorted_idx);

    % Remap paths with new sorted color indices
    cmap = lines(length(groups));
    path_color_map = containers.Map('KeyType','int32','ValueType','any');
    for i = 1:length(groups)
        for pt_id = groups{i}
            path_color_map(pt_id) = cmap(i,:);
        end
    end

    % Plot
    figure;
    hold on;
    grid on;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    scatter3(P(:,1), P(:,2), P(:,3), ...
             10*peak, 'r', 'filled');  % Highlighted markers

    
    for cid = keys(connections)
        child_id = cid{1};
        parent_ids = connections(child_id);
        child_point = P(P(:,4)==child_id, 1:3);
        for j = 1:length(parent_ids)
            parent_point = P(P(:,4)==parent_ids(j), 1:3);
            if isKey(path_color_map, child_id)
                clr = path_color_map(child_id);
                plot3([child_point(1), parent_point(1)], ...
                      [child_point(2), parent_point(2)], ...
                      [child_point(3), parent_point(3)], ...
                      '-', 'Color', clr, 'LineWidth', 2);
            end
        end
    end

    hold off;
end

% divide the groups into subgroups.
function subgroups = divideGroupsFourFold(P, groups)
    % Assumes P is [x y z idx]
    % groups: cell array of point indices

    nGroups = length(groups);
    angles = [0, pi/2, pi, 3*pi/2];  % 0°, 90°, 180°, 270°
    group_centroids = zeros(nGroups, 2);

    % Compute centroids for each group
    for i = 1:nGroups
        pts = P(groups{i}, 1:2);
        group_centroids(i,:) = mean(pts, 1);
    end

    % Center data around origin for symmetry comparison
    center = mean(group_centroids, 1);
    group_centroids_centered = group_centroids - center;

    % Rotate all centroids to get "canonical" representatives
    canonical_forms = zeros(nGroups, 2, 4);
    for a = 1:4
        R = [cos(angles(a)), -sin(angles(a));
             sin(angles(a)),  cos(angles(a))];
        canonical_forms(:,:,a) = (R * group_centroids_centered')';
    end

    % Use angle to assign groups to symmetry quadrant
    subgroups = cell(1,4);
    for i = 1:nGroups
        x = group_centroids_centered(i,1);
        y = group_centroids_centered(i,2);
        angle = mod(atan2(y, x), 2*pi);  % angle from center

        if angle < pi/2
            idx = 1;
        elseif angle < pi
            idx = 2;
        elseif angle < 3*pi/2
            idx = 3;
        else
            idx = 4;
        end

        subgroups{idx}{end+1} = groups{i}; %#ok<AGROW>
    end
end

% try to fit the local minima.
function [sigma, fitSuccess] = fitSymmetricGaussian2DAtPoint(Z, coord, window_size, nf)
% fitSymmetricGaussian2DAtPoint fits a symmetric 2D Gaussian to a local minimum in a 2D matrix.
%
% Inputs:
%   Z           - 2D matrix (e.g., intensity or elevation)
%   coord       - [x y] coordinate of the local minimum (column, row format)
%   window_size - size of square region around the point (must be odd)
%
% Outputs:
%   sigma       - estimated Gaussian standard deviation (same in x and y)
%   fitSuccess  - true if fitting was successful

    fitSuccess = false;
    sigma = NaN;

    if mod(window_size, 2) == 0
        error('window_size must be odd');
    end

    x0 = coord(1);
    y0 = coord(2);
    half_win = floor(window_size / 2);
    [rows, cols] = size(Z);

    if x0 <= half_win || x0 + half_win > cols || ...
       y0 <= half_win || y0 + half_win > rows
        warning('Window exceeds matrix bounds.');
        return;
    end

    % Extract subregion
    Z_sub = Z(y0-half_win:y0+half_win, x0-half_win:x0+half_win);
    Z_sub = nfold(Z_sub, nf);
    [X, Y] = meshgrid(1:window_size, 1:window_size);
    xdata = [X(:), Y(:)];
    zdata = normalize(Z_sub(:));

    % Initial guess: [A, x0, y0, sigma, B]
    A_init = min(zdata) - max(zdata);  % peak depth
    guess = [A_init, half_win+1, half_win+1, 0.0000001, mean(zdata)];

    % Symmetric 2D Gaussian model
    gaussSym = @(p, xdata) ...
        p(1) * exp( -((xdata(:,1)-p(2)).^2 + (xdata(:,2)-p(3)).^2) / (2*p(4)^2) ) + p(5);

    options = optimoptions('lsqcurvefit', 'Display', 'off');
    lb = [-Inf, 0, 0, 0.1, -Inf];
    ub = [0, window_size, window_size, window_size, Inf];

    try
        p_fit = lsqcurvefit(gaussSym, guess, xdata, zdata, lb, ub, options);
        sigma = p_fit(4);
        fitSuccess = true;
    catch
        warning('Fitting failed.');
    end
end

function out = nfold(in, nf)
out = in;
for i = 2:nf
    out = out + imrotate(in, (i - 1)*360/nf, 'nearest','crop');
end
out = out/nf;
end