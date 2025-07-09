%% input
% scale = 5;  % from drft_determination.m
% edge = 10;  % from drft_determination.m
% parsize = 10; % diameter, unit: nm
% 
% res_tgt = 0.5; % target xy resolution, unit: A/pix
pp2 = lpmap2_pp + 1;  % drift corrected particles [fast, slow], start from 1;
pplist = []; % particle list [raw X, raw Y, final X, final Y]
%% pre-process
% pp2(:, 2) = size(lpmap2, 1) - pp2(:, 2) + 1; % correct for the scan origin (from [0, 0] to [0, d1] for [slow, fast])

drft2_vf = drft_vf/(60*res); % unit: pix/s
drft2_vs = drft_vs/(60*res); % unit: pix/s
parsize2 = scale*parsize/res; % diameter, unit: pix (after scaling)

% [YY, XX] = ndgrid(d1:-1:1, 1:d2);
[YY, XX] = ndgrid(1:d1, 1:d2);
TT = (YY - 1)*scan_vs + (XX - 1)*scan_vf;

[YY1, XX1] = ndgrid(1:d1, 1:d2);

% drift correction
YY2 = scale*((YY-1) - TT*drft2_vs);   % start from 0
XX2 = scale*((XX-1) - TT*drft2_vf);   % start from 0

YY3 = YY2 + scale*edge + 1;
XX3 = XX2 + scale*edge + 1;

%% optional - locate picked particles from the drift corrected map
%%% this is where the particles should be picked
%%% this is only for display purpose
[YYY, XXX] = ndgrid(1:size(lpmap2, 1), 1:size(lpmap2, 2));
canvas = YYY.*0;
canvas_ct = canvas.*0;
for i = 1:size(pp2, 1)
    dd = sqrt((XXX - pp2(i, 1)).^2 + (YYY - pp2(i, 2)).^2);
    sel = dd <= parsize2/2;
    canvas = canvas + lpmap2.*sel;
    canvas_ct = canvas_ct + sel;
    % MIJ.createImage(lpmap2.*sel)
end
canvas = canvas./canvas_ct;
% canvas(isnan(canvas)) = min(canvas(:));
MIJ.createImage(canvas);

%% optional - locate picked particles from the raw map
%%% this is only for display purpose
canvas = zeros(size(lpmap));
canvas_ct = canvas.*0;
for i = 1:size(pp2, 1)
    dd = sqrt((XX3 - pp2(i, 1)).^2 + (YY3 - pp2(i, 2)).^2);
    sel = dd <= parsize2/2;

    lpmap_sel = lpmap(sel);
    canvas = canvas + lpmap.*sel;
    canvas_ct = canvas_ct + sel;
    % MIJ.createImage(lpmap.*sel)
end
canvas = canvas./canvas_ct;
% canvas(isnan(canvas)) = min(canvas(:));
MIJ.createImage(canvas);

%% optional - locate picked particles from the raw map to the new space (drift corrected space)
%%% this is only for display purpose
canvas = zeros(size(lpmap2));
canvas_ct = canvas.*0;
for i = 1:size(pp2, 1)
    dd = sqrt((XX3 - pp2(i, 1)).^2 + (YY3 - pp2(i, 2)).^2);
    sel = dd <= parsize2/2;

    YY3_sel = YY3(sel);
    XX3_sel = XX3(sel);
    YY1_sel = YY1(sel);
    XX1_sel = XX1(sel);

    lpmap_sel = lpmap(sel);
   
    canvasi = canvas.*0;
    for j = 1:numel(YY3_sel)
         % canvasi(round(YY3_sel(j)), round(XX3_sel(j))) = lpmap_sel(j);
         canvasi(round(YY3_sel(j)), round(XX3_sel(j))) = lpmap(YY1_sel(j), XX1_sel(j));
         canvas_ct(round(YY3_sel(j)), round(XX3_sel(j))) = 1;
    end
    canvas = canvas + canvasi;
    % MIJ.createImage(canvasi)
end
canvas = canvas./canvas_ct;
% canvas(isnan(canvas)) = min(canvas(:));
MIJ.createImage(canvas);   % note: this is vertically transformed for display!

%% output drift corrected data aquasition locations for single particls
%%% this is drift corrected unaligned x,y coordiates of the phase data aquasition
%%% locations for individual particles, output for 3D construction
scale2 = res*10/res_tgt; 
YY2 = scale2*((YY-1) - TT*drft2_vs);   % start from 0
XX2 = scale2*((XX-1) - TT*drft2_vf);   % start from 0

pplist = [];
idx = numel(pplist);
ppmap = zeros(floor(10*parsize/res_tgt), floor(10*parsize/res_tgt), size(pp2, 1));
ppmap_ct = ppmap.*0;
for i = 1:size(pp2, 1)
    idx = idx + 1;
    dd = sqrt((XX3 - pp2(i, 1)).^2 + (YY3 - pp2(i, 2)).^2);
    sel = dd < 0.99*parsize2/2;
    YY2_sel = YY2(sel);
    XX2_sel = XX2(sel);
    YY1_sel = YY1(sel);
    XX1_sel = XX1(sel);


    YY2_norm = YY2 - (scale2/scale)*(pp2(i, 2) - 1 - scale*edge);
    XX2_norm = XX2 - (scale2/scale)*(pp2(i, 1) - 1 - scale*edge);

    YY2_sel_norm = YY2_norm(sel);
    XX2_sel_norm = XX2_norm(sel);

    YY2_sel_disp = YY2_sel_norm + (5*parsize/res_tgt);
    XX2_sel_disp = XX2_sel_norm + (5*parsize/res_tgt);

    for j = 1:numel(YY2_sel_disp)
        ppmap(round(YY2_sel_disp(j) + 1), round(XX2_sel_disp(j) + 1), i) = lpmap(YY1_sel(j), XX1_sel(j));
        ppmap_ct(round(YY2_sel_disp(j) + 1), round(XX2_sel_disp(j) + 1), i) = 1;
    end

    
    YY3 = YY2_norm*res_tgt/10;   % unit: nm
    XX3 = XX2_norm*res_tgt/10;   % unit: nm
    pplist{idx} = [YY1(:) XX1(:) YY3(:) XX3(:)]; 
end

ppmap = ppmap./ppmap_ct;
ppmap(isnan(ppmap)) = min(ppmap(:));
MIJ.createImage(ppmap);

%% write output file
% save(filename + "_pplist.mat", "edge", "scale", "lpmap2", "lpmap2_pp", "parsize", "pplist", "res_tgt". "-append");
