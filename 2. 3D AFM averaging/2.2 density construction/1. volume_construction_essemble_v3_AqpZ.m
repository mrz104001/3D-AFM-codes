%% input
filename = "structset_AqpZ_test5_0807a";
% filename = "structset_AqpZ_test4_0912a";
% filename = "structset_AqpZ_0722a";
t = datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ssZ');
filesuffix = string(t);
parsize = 10; % unit: nm

xyres_tgt = 0.4;  % unit: A/pix;
% zres_tgt = xyres_tgt;   % unit: A/pix
zres_tgt = 0.065;   % unit: A/pix;
nf = 4;
% xysig = 2.75;   % unit: A
% zsig = 0.3;    % unit: A
xysig = 2.75;   % unit: A
zsig = 0.75;    % unit: A
%% particle selection
load(filename + ".mat")
parlist = zeros(numel(structset)*10, 7);
idx = 0;
for s = 1:numel(structset)
    ss = structset{s};
    palist = ss.palist;
    for p = 1:size(palist, 1)
        idx = idx + 1;
        parlist(idx, 1) = s;
        parlist(idx, 2) = p;
        parlist(idx, 3:end) = palist(p, :);
    end

end
parlist = parlist(1:idx, :);

%% volume construction
%%%%
threshccv = 0.9;
threshisv = 0.6;
% threshccv = 0;
% threshisv = 0;
ss_tbp = 1:numel(structset);
% ss_tbp = 3:4;

%%%%
htmax = -Inf;
htmin = Inf;
for s = 1:numel(structset)
    ss = structset{s};
    ht3 = ss.matHT;
    pplist = ss.pplist;
    sel = ht3(:, :, 1).*0;
    for i = 1:numel(pplist)
        pp = pplist{i};
        for j = 1:size(pp, 1)
            sel(pp(j, 1), pp(j, 2)) = 1;
        end
    end
    sel = sel > 0;
    bg = mean(ht3(:, :, 1:3), 3);
    ht4 = ht3 - prctile(bg(sel), 50);   % background adjustment
    htmax = max([max(ht4(:)), htmax]);
    htmin = min([min(ht4(:)), htmin]);
end



ppvolume = zeros(floor(10*parsize/xyres_tgt), floor(10*parsize/xyres_tgt), floor(10*(htmax-htmin)/zres_tgt)+1);
ppvolume_ct = ppvolume.*0;

ct_tot = 0;
ct_bad = 0;

for s = 1:numel(structset)
    if ismember(s, ss_tbp)
        ss = structset{s};
        pplist = ss.pplist2;
        palist = ss.palist;
        ht3 = ss.matHT;
        ee4 = ss.matEN;
        % ee4 = max(ee4(:)) - ee4;  % 240807, for direct calculation of water density
        pr = exp(-ee4);
        sel = ht3(:, :, 1).*0;
        for i = 1:numel(pplist)
            pp = pplist{i};
            for j = 1:size(pp, 1)
                sel(pp(j, 1), pp(j, 2)) = 1;
            end
        end
        sel = sel > 0;
        bg = mean(ht3(:, :, 1:3), 3);
        ht4 = ht3 - prctile(bg(sel), 50);   % background adjustment
        ht4 = 10*(ht4 - htmin)./zres_tgt;  % convert from nm to pix (res_tgt has a unit of A/pix)
        for i = 1:numel(pplist)
            ppi = pplist{i};
            if palist(i, 4) > threshccv && palist(i, 5) > threshisv

                YY1_sel = ppi(:, 1);
                XX1_sel = ppi(:, 2);
                YY2_sel_norm = 10*ppi(:, 3)/xyres_tgt;
                XX2_sel_norm = 10*ppi(:, 4)/xyres_tgt;
                YY2_sel_disp = YY2_sel_norm + (5*parsize/xyres_tgt);
                XX2_sel_disp = XX2_sel_norm + (5*parsize/xyres_tgt);

                for j = 1:size(ppi, 1)
                    flag = sqrt(ppi(j, 3).^2 + ppi(j, 4).^2) < parsize/2;
                    for k = 1:size(pr, 3)  
                        if ~isnan(pr(YY1_sel(j), XX1_sel(j), k)) && flag
                            ct_tot = ct_tot + 1;
                            try
                                ppvolume(round(YY2_sel_disp(j) + 1), round(XX2_sel_disp(j) + 1),...
                                    round(ht4(YY1_sel(j), XX1_sel(j), k))+1) = ...
                                    ppvolume(round(YY2_sel_disp(j) + 1), round(XX2_sel_disp(j) + 1),...
                                    round(ht4(YY1_sel(j), XX1_sel(j), k))+1) + ...
                                    pr(YY1_sel(j), XX1_sel(j), k);
                                ppvolume_ct(round(YY2_sel_disp(j) + 1), round(XX2_sel_disp(j) + 1),...
                                    round(ht4(YY1_sel(j), XX1_sel(j), k))+1) = ...
                                    ppvolume_ct(round(YY2_sel_disp(j) + 1), round(XX2_sel_disp(j) + 1),...
                                    round(ht4(YY1_sel(j), XX1_sel(j), k))+1) + 1;
                            catch exception
                                ct_bad = ct_bad + 1;
                                % warning("data point out of the boundary");
                            end
                        end
                    end
                end
            end
        end
    end
end

disp("volume constructed... " + string(100 - 100*ct_bad/ct_tot) + "% data point used...");
ppvolume_ct2 = ppvolume_ct;
ppvolume_ct(ppvolume_ct == 0) = 1;
ppvolume = ppvolume./ppvolume_ct;
%%
%%%% select detections
trust_zz = 1100;   % unit: pixel
[XX, YY, ZZ] = ndgrid(1:size(ppvolume, 1), 1:size(ppvolume, 2), 1:size(ppvolume, 3));
dd = sqrt((XX - size(ppvolume, 1)/2).^2 + (YY - size(ppvolume, 2)/2).^2);
sel = dd < size(ppvolume, 1)/2 - 10;
% sel = dd < size(ppvolume, 1)/2 - 100;
selz =  ZZ < trust_zz;
%%%% construct detection volume
ppvolume2 = ppvolume.*sel.*selz;
ppvolume2 = nfold(ppvolume2, nf);
%%%% construct density volume
%%
sigxy = xysig/xyres_tgt;   % unit: pix
sigz = zsig/zres_tgt;   % unit: pix
h = make_3D_LAFM_kernel1a(sigxy, sigz);
% h = make_3D_LAFM_kernel1a(5, 2);
ppvolume3 = imfilter(ppvolume2, h);
ppvolume3 = nfold(ppvolume3, nf);
ppvolume3(~sel) = min(ppvolume3(:));
%%%% display
% MIJ.createImage(ppvolume2);
try
    MIJ.createImage(ppvolume3);
catch
    disp('No MIJ plugin detected, optional plotting skipped.')
end

%% write output files

save(filename + "_v3_" + filesuffix + ".mat", "ppvolume2", "ppvolume3", "zsig", "xysig", "zres_tgt", "xyres_tgt", "nf", "parsize");
disp("file: " + filename + "_v3_" + filesuffix + ".mat created...");
%% helper function
function output = nfold(input, nf)
output = input;
for i = 2:nf
    output = output + imrotate3(input, (i-1)*(360/nf), [0 0 1]);
end
output = output/nf;
end
