%% find particle coordinate from the lowest point mapping data (lpm)
%%%% output lpm data for as the input data for Particle Picker
lpmap = -1*lpm;
[d1, d2, d3] = size(lpmap);
MIJ.createImage(lpmap);
filename = "AqpZ 2D scan test 2 trace";
energydata = E;
heightdata = z_avg_m_corr_filled_sub;
phasedata = phase_avg_m_corr;
forcedata = F;
ampdata = r_avg_m_corr;
%%%% write output file
save(filename + "_input.mat", "lpmap", "filename", "energydata", "heightdata", "d1", "d2", "d3", "phasedata", "forcedata", "ampdata");
disp("new file: " + filename + "_input.mat" + " created...")
%% input coordinates from Particle Picker (ImageJ)
lpmap_pp;  % read particle coordinates
lpmap_np;  % read neighbor pair information
% lpmap_pp = [19 20; 19 22; 22 19;24 21; 21 24; 24 24];
% lpmap_np = [0 1; 0 2; 2 3; 1 4; 4 5; 3 5];

%%%% write output file
save(filename + "_input.mat", "lpmap_pp", "lpmap_np", "-append");
disp("existed file: " + filename + "_input.mat" + " updated...")
disp("ready for drift correction...")
%% drift correction
drft_vf0 = -3; % unit: nm/min
drft_vs0 = 1; % unit: nm/min
scan_vf = 0.025; % unit: s/pix;
scan_vs = 3; % unit: s/pix;
res = 0.5; % unit: nm/pix
rng(1);

%%%% run script
run("drift_determination_v1.m")
disp("drift correction completed...")
%%%% write output file
save(filename + "_drft.mat", "drft_vf", "drft_vf0", "drft_vs", "drft_vs0", "res", "scan_vf", "scan_vs", "lpmap2");
disp("new file: " + filename + "_drft.mat" + " created...")
disp("next: input coordinates for particle determination...")
%% input coordinates for particle determination from Particle Picker (ImageJ)
%%%% input
lpmap2_pp;  % read drift corrected scaled particle coordinates from Particle Picker
% lpmap2_pp = [247 128; 239 218; 335 136; 329 226];
%%%% write output file
save(filename + "_pplist.mat", "lpmap2", "lpmap2_pp");
disp("new file: " + filename + "_pplist.mat" + " created...")
disp("ready for particle determination...")
%% particle determination
%%%% input
parsize = 10; % diameter, unit: nm
res_tgt = 0.5; % target xy resolution, unit: A/pix

%%%% run script
% run("particle_determination_v1.m")    % this script output coordinates in unit: pixel at res_tgt
% run("particle_determination_v2.m")    % this script output coordinates in unit: nm
run("particle_determination_v2b.m")    % this script output ALL coordinates in unit: nm
disp("particle determination completed...")
%%%% write output file
save(filename + "_pplist.mat", "edge", "scale", "parsize", "pplist", "res_tgt", "-append");
disp("existed file: " + filename + "_pplist.mat" + " updated...")
disp("unaglined particle x,y coordinates of data: " + filename +  ". is determined...")

%% energy background processing
%%%% this is not built no previous resutls
%%%% input
zconst = -12.6*2; % from nm/V
ee = energydata;  % corrected 0719 8:30pm
ht = heightdata;


%%%% run script
% run("energy_background_v1.m");
% run("energy_background_v2.m");
run("energy_background_v3.m");   % starting 080624, remove background alignment
disp("energy background processing completed...")
%%%% write output file
save(filename + "_eeht.mat", "zconst", "ee4", "ee3", "ht3");
disp("new file: " + filename + "_eeht.mat" + " created...")
disp("energy background of data: " + filename + ". is ready for 3D volume construction...")

%% update structset data
%%%% read input file
stuctsetname = "structset_AqpZ_test7_0828a";
load(stuctsetname + ".mat");

%%%%
s = struct;
s.filename = filename;
s.mapLP = lpmap2;
s.pplocs = lpmap2_pp;
s.pplist = pplist;
s.matEN = ee4;
s.matHT = ht3;

idx = numel(structset) + 1;
% idx = 1;
structset{idx} = s;
disp("structset updated...")
%%%% write output file
save(stuctsetname + ".mat", "structset", "-append");
disp("file: " + stuctsetname + ".mat" + " updated...")