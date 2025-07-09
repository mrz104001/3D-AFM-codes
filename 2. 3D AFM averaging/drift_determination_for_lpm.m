%% input
% manual input of the particle position and relations.
% need four points in the square lattice 2D patch.
clear
close all
clc
lpmap_pp=[17,16;17,38;35,19;36,40];

frame_time=120;   % unit: s
pixels=60;      % unit;
size_nm=30;     % unit: nm. scan frame size.

scan_vf = frame_time/2/pixels^2; % unit: s/pix;
scan_vs = frame_time/pixels; % unit: s/pix;
res = size_nm/pixels; % unit: nm/pix

% need to flip the coordinates in y, the coordinate in ImageJ is inverted
% in y direction.
lpmap_pp(:,2)=pixels-1-lpmap_pp(:,2);
%% pre-process

% correct for the scan origin (from [0, 0] to [0, d1] for [slow, fast])
% pp(:, 2) = d1 - pp(:, 2) + 1;

% scan elipse time
tt =lpmap_pp(:, 1)*scan_vf + lpmap_pp(:, 2)*scan_vs;
%% drift speed search
close all
drft_allow_s = res/scan_vs*10;  % unit: nm/s, half of the slow axis scan speed.


drft_range_s = -drft_allow_s:drft_allow_s/100:drft_allow_s;


drft_range_s = drft_range_s./(res);


[I, J] = ndgrid(drft_range_s, drft_range_s);
VV = I.*0;
VV_r=VV;
VV_phase=VV;
pp_ij = lpmap_pp;
% warning('off')
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        pp_ij = pp_ij.*0;
        pp_ij(:, 1) = lpmap_pp(:, 1) + tt*I(i, j);
        pp_ij(:, 2) = lpmap_pp(:, 2) + tt*J(i, j);
        
        % calculating the polar coordinates with respect to the center.
        pp_ij(:,1)=pp_ij(:,1)-mean(pp_ij(:,1),"all");
        pp_ij(:,2)=pp_ij(:,2)-mean(pp_ij(:,2),"all");

        [pp_ij_phase,pp_ij_r]=cart2pol(pp_ij(:,1),pp_ij(:,2));
        pp_ij_phase=pp_ij_phase/pi*180;

        % sort the matrix to get the neighboring ones.
        if i==1&&j==1
            [~,ind]=sort(pp_ij_phase);
        end
        VV_r(i,j)=abs(pp_ij_r(ind(1))-pp_ij_r(ind(2)));
        VV_phase(i,j)=abs(pp_ij_phase(ind(1))-pp_ij_phase(ind(2)));
        
        %%% ellipise fitting strategy

    end
end

VV_phase=abs(VV_phase-90);

VV_r=abs(VV_r);

VV=VV_r./max(VV_r,[],'all')+VV_phase./max(VV_phase,[],'all');
figure;
imagesc([-drft_allow_s drft_allow_s],[-drft_allow_s drft_allow_s],VV);
shading flat
axis tight square xy
title('Square point for xy drift determination');
ylabel('fast axis drift speed(x, nm/s)');
xlabel('slow axis drift speed(y, nm/s)');

%% test for two sets of different search range.

close all
drft_allow_s = 1/scan_vs/2;  % unit: pixels/s, half of the slow axis scan speed.
drft_allow_f = 1/scan_vf/pixels;  % unit: pixels/s, half of the fast axis scan speed.

drft_range_s = -drft_allow_s:drft_allow_s/100:drft_allow_s;
drft_range_f = -drft_allow_f:drft_allow_f/100:drft_allow_f;

[I, J] = ndgrid(drft_range_f, drft_range_s);
VV = I.*0;
VV_r=VV;
VV_phase=VV;
pp_ij = lpmap_pp;
% warning('off')
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        pp_ij = pp_ij.*0;
        pp_ij(:, 1) = lpmap_pp(:, 1) + tt*I(i, j);
        pp_ij(:, 2) = lpmap_pp(:, 2) + tt*J(i, j);
        
        % calculating the polar coordinates with respect to the center.
        pp_ij(:,1)=pp_ij(:,1)-mean(pp_ij(:,1),"all");
        pp_ij(:,2)=pp_ij(:,2)-mean(pp_ij(:,2),"all");

        [pp_ij_phase,pp_ij_r]=cart2pol(pp_ij(:,1),pp_ij(:,2));
        pp_ij_phase=pp_ij_phase/pi*180;

        % sort the matrix to get the neighboring ones.
        if i==1&&j==1
            [~,ind]=sort(pp_ij_phase);
        end
        VV_r(i,j)=abs(pp_ij_r(ind(1))-pp_ij_r(ind(2)));
        VV_phase(i,j)=abs(pp_ij_phase(ind(1))-pp_ij_phase(ind(2)));
        
        %%% ellipise fitting strategy

    end
end

VV_phase=abs(VV_phase-90);

VV_r=abs(VV_r);

VV=VV_r./max(VV_r,[],'all')+VV_phase./max(VV_phase,[],'all');
figure;
imagesc([-drft_allow_s drft_allow_s].*res,[-drft_allow_f drft_allow_f].*res,VV);
shading flat
axis tight square xy
title('Square point for xy drift determination');
ylabel('fast axis drift speed(x, nm/s)');
xlabel('slow axis drift speed(y, nm/s)');

%% drift speed
[i, j] = find(VV == min(VV(:)));
drft2_vf =I(i, j); % unit: pix/s
drft2_vs =J(i, j); % unit: pix/s
drft_vf = drft2_vf*res; % unit: nm/s
drft_vs = drft2_vs*res; % unit: nm/s

if (drft_vf<1||drft_vs<1)
    disp("fast-(x-)axis drift speed " + drft_vf*60 + " nm/min")
    disp("slow-(y-)axis drift speed " + drft_vs*60 + " nm/min");
else
    disp("fast-(x-)axis drift speed " + drft_vf + " nm/s")
    disp("slow-(y-)axis drift speed " + drft_vs + " nm/s");
end

%% optional - display particle locations
%%% this is only for display purpose

% pp_ij(:, 1) = pp(:, 1) - tt*(drft2_vf0+I(i, j));
% pp_ij(:, 2) = pp(:, 2) - tt*(drft2_vs0+J(i, j));
% dd_ij = pp_ij(nn(:, 1), :) - pp_ij(nn(:, 2), :);
% dd_ij = vecnorm(dd_ij')';
% canvas = lpmap.*0;
% pp2 = pp;
% pp2 = round(pp_ij);
% % lpmap_pp2 = pp2;
% % lpmap_pp2(:, 2) = d1 - lpmap_pp2(:, 2) + 1;
% for k = 1:size(pp, 1)
%     canvas(pp2(k, 2) + 1, pp2(k, 1) + 1) = 1;
% end
% MIJ.createImage(canvas)
% %% output drift corrected data
% %%% this is the map where particles should be picked for further analysis
% % [YY, XX] = ndgrid(d1:-1:1, 1:d2);
% [YY, XX] = ndgrid(1:d1, 1:d2);
% TT = (YY - 1)*scan_vs + (XX - 1)*scan_vf;
% 
% scale = 5;
% edge = 10;
% YY2 = round(scale*((YY-1) - TT*drft2_vs));  % start from 0
% XX2 = round(scale*((XX-1) - TT*drft2_vf));  % start from 0
% 
% canvas = zeros(scale*(d1+2*edge), scale*(d2+2*edge));
% data = lpmap - min(lpmap(:));
% for i = 1:d1
%     for j = 1:d2
%         canvas(scale*edge + YY2(i, j) + 1, scale*edge + XX2(i, j) + 1) = data(i, j);
%     end
% end
% canvas(canvas == 0) = nan;
% canvas = fillmissing2(canvas,'linear');
% MIJ.createImage(canvas)  % note: this is vertically transformed for display!
% 
% lpmap2 = canvas;
%% give corrected image.
% drift_x=drft_vf;           % input the drift speed of x and y, nm/s.
% drift_y=drft_vs;

drift_x=0.025;           % input the drift speed of x and y, nm/s.
drift_y=-0.45;

if (~exist('lpm','var'))          % ask for new data if the first section is called to clear all data in workspace.
    % must define the format of the file to be read.
    file_format='.tif';

    % read the file based on the previous loading and file format.
    [name0,path0]=uigetfile({['*',file_format]},['Lowest point mapping, choose ',file_format,' file to read'],path);
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

    lpm=tiffreadVolume(path);
    lpm=flip(lpm,1);
end

% use the above information to calculate the real coordinate for the data
% points.


    % first, give a time for each points.
    time_delay=[];      % record the time needed to get to the current point from the start, in s.
    for i=1:size(lpm,1)     % y direction 
        for j=1:size(lpm,2)         % x direction
            time_delay(i,j)=(j - 1)*scan_vf + (i - 1)*scan_vs;
        end
    end

    % to calculate the real coordinate offset for each point.
    offset_x=time_delay*drift_x;
    offset_y=time_delay*drift_y;

    [o_x,o_y]=meshgrid(1:size(lpm,2),1:size(lpm,1));
    o_x=o_x*res;
    o_y=o_y*res;

    real_x=o_x+offset_x;
    real_y=o_y+offset_y;


    % try to plot the data.
    N=pixels;
    xlin=linspace(min(real_x(:)),max(real_x(:)),N);
    ylin=linspace(min(real_y(:)),max(real_y(:)),N);
    [X,Y]=meshgrid(xlin,ylin);
    Z=griddata(real_x(:),real_y(:),double(lpm(:)),X,Y,'natural');
    % figure;
    % imagesc(Z);
    % axis equal tight;


MIJ.createImage(flip(Z,1));
%% write output file
% save(filename + "_drft.mat", "drft_vf", "drft_vf0", "drft_vs", "drft_vs0", "res", "scan_vf", "scan_vs". "lpmap2");