%% input
% zconst = -12.6*2; % from nm/V
% ee = energydata;
% ht = heightdata;
%% height data process - construction
% background correction (z-drift correction)
ht2 = ht;
% unit conversion (from V to nm)
ht2 = ht2 * zconst;
% match to the energy data
ht3 = ht2(:, :, 1:size(ee, 3));

%% phase data background 2 - construction
econst = 1/4.11E-21; % unit: J/kBT
UpLim = [10 1000 10 10 1000];
LoLim = [0 0 -10 0 0];
StPoint = [1 1 0 1 1];
foption = fitoptions('method','NonlinearLeastSquares','Upper',UpLim,...
        'Lower', LoLim, 'StartPoint', StPoint, 'normalize', 'off');
fittype1 = fittype('a * exp(-x / b) + d * exp(-x / e) + c ',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a' 'b' 'c' 'd' 'e'});

ee2 = ee - min(ee(:));
ee3 = ee2.*0;
coeffs = cell(d1, d2);

xfs = cell(d1, d2);
yfs = cell(d1, d2);
exps0s = cell(d1, d2);


for i = 1:d1
    disp("progress: " + string(100*i/d1) + "%...");
    for j = 1:d2
        % ee2(i, j, :) = flip(ee2(i, j, :));
        xf = squeeze(ht3(i, j, :));
        yf = squeeze(ee2(i, j, :));
        [yf, nc, ns] = normalize(yf);
        try
        exps0 = fit(xf(~isnan(yf)), yf(~isnan(yf)),fittype1,foption);
        coeffs0 = coeffvalues(exps0);
        coeffs{i, j} = coeffs0;
        ee3n = yf - (coeffs0(1).*exp(-xf./coeffs0(2)) + coeffs0(4).*exp(-xf./coeffs0(5)) + coeffs0(3));
        % ee3(i, j, :) = ee3n .* ns + nc;
        ee3(i, j, :) = ee3n .* ns;
        catch exception
            warning("Error found...");
            ee3(i, j, :) = nan;
        end

        xfs{i, j} = xf;
        yfs{i, j} = yf;
        exps0s{i, j} = exps0;
    end
end
% ee3 = ee2;  % delete
ee4 = ee3 * econst;

%% write output files
% save(filename + "_eeht.mat", "zconst", "ee4", "ht3");
%% optional - display background fitting and bg-removed phase data
%%% this is only for display purpose

sel = zeros(d1, d2, d3);
for i = 1:numel(pplist)
    ppi = pplist{i};
    for j = 1:size(ppi, 1)
        sel(ppi(j, 1), ppi(j, 2), :) = 1;
    end
end
% sel = sel == 0;


xf = xfs{10, 10};
yf = yfs{10, 10};
exps0 = exps0s{10, 10};
figure()
hold on
plot(xf(~isnan(yf)), yf(~isnan(yf)), "o");
plot(exps0)
% ph3 = ph;
figure()
hold on
for i = 1:d1
    for j = 1:d2
        if sel(i, j, 1)
            % plot(squeeze(ht3(i,j,:)), squeeze(ee4(i,j,:)));
            plot(squeeze(ee4(i,j,:)));
        end
    end
end
hold off
