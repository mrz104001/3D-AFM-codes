filename = "structset_AqpZ_0722a_v3_8-Jul-2025-17-14-43-0400";
load(filename + ".mat");
in_ppvolume = ppvolume4;
in_ppvolume(in_ppvolume == 0) = min(in_ppvolume(in_ppvolume>0));
%%
econst = 4.11E-21; % unit: kBT/J
% ppvolumeE = log(in_ppvolume);   % energy, unit: kBT, for direct calculation of force from water density
ppvolumeE = -log(in_ppvolume);   % energy, unit: kBT
ppvolumeF = -(ppvolumeE(:, :, 2:end) - ppvolumeE(:, :, 1:end-1)); 
% ppvolumeF = ppvolumeF./0.02E-9;   % force, unit: kBT/m
ppvolumeF = ppvolumeF./(zres_tgt*1E-10);   % force, unit: kBT/m   %082724
ppvolumeF = ppvolumeF.*econst; % force, unit: N
%%
save(filename + ".mat", "in_ppvolume", "ppvolumeE", "ppvolumeF", "-append"); 
disp("file:" + filename + ".mat updated...");
% 

%%
% filename = "structset_AqpZ_test5_0807a_v3_F2_14-Aug-2024-14-07-48-0400";
% load(filename + ".mat");


F2_z_base_A = 1.3;  % unit: angstrom  % AqpZ sm, AqpZ lattice, mica 
F2_z_base_A = 4.3;

% F2_z_len_A = 13.8;  % unit: angstrom   % mica lattice
F2_z_len_A = 19.8;  % unit: angstrom   % Aqpz sm, AqpZ lattice

% F2_z_base = 181;   % unit: voxel
% F2_z_top = 660;   % unit: voxel
F2_z_base = round(F2_z_base_A/zres_tgt) + 1;   % unit: voxel
F2_z_top = round((F2_z_base_A + F2_z_len_A)/zres_tgt);   % unit: voxel


ppvolumeF2 = ppvolumeF(:, :, F2_z_base:F2_z_top);
sz = size(ppvolumeF2);
[~, ~, ZZ] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));

ZZ2 = ZZ.*0;

for i = 1:sz(1)
    for j = 1:sz(2)
        ppcol = squeeze(ppvolumeF2(i, j, :));
        ppmax = find(islocalmax(ppcol), 1, "first");
        if isempty(ppmax)
            ppmax = sz(3);
        end
        ZZ2(i, j, :) = ppmax;
    end
end 


sel = ZZ < ZZ2;
ppvolumeF2(sel) = nan;

%%
save(filename + ".mat", "ppvolumeF2", "F2_z_base", "-append");
disp("file:" + filename + ".mat updated...");

