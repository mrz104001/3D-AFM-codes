
sigxy = xysig/xyres_tgt;   % unit: pix
sigz = zsig/zres_tgt;   % unit: pix
h = make_3D_LAFM_kernel1a(sigxy, sigz);

% sz = size(ppvolumeF2);
% h1 = ones(sz);
% h2 = imfilter(h1, h);
% 
% sc = sum(h1(:))/sum(h2(:));

sz = size(h);
md = 0.5*(sz + 1);
rad = sz/2;
rad = rad*2.355/3;   % FWHM equivalent radius
sc = 4*pi*rad(1)*rad(2)*rad(3)*h(md(1),md(2),md(3))/3;

ppvolumeF3 = ppvolumeF2*sc;


%%
save(filename + ".mat", "sc", "ppvolumeF3", "-append");
disp("file:" + filename + ".mat updated...");