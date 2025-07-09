%% input
% manual input of the particle position and relations.
drft_vf0 = -1; % unit: nm/min
drft_vs0 = 1; % unit: nm/min

frame_time=1;   % unit: s
pixels=60;      % unit;

scan_vf = 0.025; % unit: s/pix;
scan_vs = 3; % unit: s/pix;
res = 0.5; % unit: nm/pix
% rng(1);
pp = lpmap_pp + 1;  % particles [fast, slow], start from 1
nn = lpmap_np + 1;  % neighbor particle pairs [particle1 particle2]
%% pre-process
% change all unit to pix/s
drft2_vf0 = drft_vf0/(60*res); % unit: pix/s
drft2_vs0 = drft_vs0/(60*res); % unit: pix/s

% correct for the scan origin (from [0, 0] to [0, d1] for [slow, fast])
% pp(:, 2) = d1 - pp(:, 2) + 1;

% scan elipse time
tt = (pp(:, 1) - 1)*scan_vf + (pp(:, 2) - 1)*scan_vs;
%% drift speed search
drft_allow = 1;  % unit: nm/min
drft_range = -drft_allow:0.01:drft_allow;
% drft_allow = 0.1;  % unit: nm/min
% drft_range = -drft_allow:0.001:drft_allow;
drft_range = drft_range./(60*res);

[I, J] = ndgrid(drft_range, drft_range);
VV = I.*0;
pp_ij = pp;
% warning('off')
for i = 1:size(I, 1)
    for j = 1:size(I, 1)
        pp_ij = pp_ij.*0;
        pp_ij(:, 1) = pp(:, 1) - tt*(drft2_vf0+I(i, j));
        pp_ij(:, 2) = pp(:, 2) - tt*(drft2_vs0+J(i, j));

        %%% ellipise fitting strategy
        pp_ij_pair0 = [pp_ij(nn(:, 1), :) - pp_ij(nn(:, 2), :);
            pp_ij(nn(:, 2), :) - pp_ij(nn(:, 1), :)];
        pp_ij_pair_mir = pp_ij_pair0;
        pp_ij_pair_mir(:, 1) = -pp_ij_pair_mir(:, 1);
        pp_ij_pair = [pp_ij_pair0; pp_ij_pair_mir];
        
        try
        ellipse_t = fit_ellipse(pp_ij_pair(:, 1),pp_ij_pair(:, 2), 0);
        ellispe_a = max([ellipse_t.a ellipse_t.b]);
        ellispe_b = min([ellipse_t.a ellipse_t.b]);
        VV(i, j) = sqrt(1 - ellispe_b^2/ellispe_a^2);
        catch expection
            VV(i, j) = nan;
        end
        %%% ellipise fitting strategy

    end
end
% warning('on')
%% drift speed
[i, j] = find(VV == min(VV(:)));
drft2_vf = drft2_vf0+I(i, j); % unit: pix/s
drft2_vs = drft2_vs0+J(i, j); % unit: pix/s
drft_vf = drft2_vf*(60*res); % unit: nm/min
drft_vs = drft2_vs*(60*res); % unit: nm/min

disp("fast-(x-)axis drift speed " + drft_vf + " nm/min")
disp("slow-(y-)axis drift speed " + drft_vs + " nm/min");

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
%% output drift corrected data
%%% this is the map where particles should be picked for further analysis
% [YY, XX] = ndgrid(d1:-1:1, 1:d2);
[YY, XX] = ndgrid(1:d1, 1:d2);
TT = (YY - 1)*scan_vs + (XX - 1)*scan_vf;

scale = 5;
edge = 10;
YY2 = round(scale*((YY-1) - TT*drft2_vs));  % start from 0
XX2 = round(scale*((XX-1) - TT*drft2_vf));  % start from 0

canvas = zeros(scale*(d1+2*edge), scale*(d2+2*edge));
data = lpmap - min(lpmap(:));
for i = 1:d1
    for j = 1:d2
        canvas(scale*edge + YY2(i, j) + 1, scale*edge + XX2(i, j) + 1) = data(i, j);
    end
end
canvas(canvas == 0) = nan;
canvas = fillmissing2(canvas,'linear');
MIJ.createImage(canvas)  % note: this is vertically transformed for display!

lpmap2 = canvas;
%% write output file
% save(filename + "_drft.mat", "drft_vf", "drft_vf0", "drft_vs", "drft_vs0", "res", "scan_vf", "scan_vs". "lpmap2");