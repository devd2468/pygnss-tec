function K = nplotKeogram(dayRoots, stations, outDir, sysList, latBand, lonSlices)
% NPLOTKEOGRAM  ROTI keograms in the Suraina/Abadi style:
%   top: zonal cross-section (time vs longitude) at an equatorial band
%   bottom: meridional slices (time vs latitude) at chosen longitudes.
% 10-min time bins; 0.5-deg spatial bins; cell value = max ROTI.
% Saves nplotKeogram.png + keogram_data.mat (reused by nplotDrift).
% Defaults: sysList={'G'}, latBand=[5 25], lonSlices=[78 88].
    if nargin < 4 || isempty(sysList), sysList = {'G'}; end
    if nargin < 5 || isempty(latBand), latBand = [5 25]; end
    if nargin < 6 || isempty(lonSlices), lonSlices = [78 88]; end
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % ---- pool IPPs across days/stations/systems ----
    AL = []; LO = []; TH = []; RO = [];  % lat, lon, absHour, roti
    for d = 1:numel(dayRoots)
        for s = 1:numel(stations)
            for t = 1:numel(sysList)
                try, P = getIPP(dayRoots{d}, stations{s}, sysList{t});
                catch, continue; end
                if isempty(P) || ~isfield(P, 'lat'), continue; end
                R = P.roti;
                for p = 1:size(R, 2)
                    v = R(:, p);
                    ok = isfinite(v) & isfinite(P.lat(:, p)) & isfinite(P.lon(:, p));
                    if sum(ok) < 5, continue; end
                    rows = find(ok);
                    AL = [AL; P.lat(rows, p)]; %#ok<AGROW>
                    LO = [LO; P.lon(rows, p)]; %#ok<AGROW>
                    TH = [TH; (d-1)*24 + (rows-1)/3600]; %#ok<AGROW>
                    RO = [RO; v(ok)]; %#ok<AGROW>
                end
            end
        end
    end
    if isempty(RO)
        fprintf('nplotKeogram: no IPP ROTI pooled.\n');
        K = [];
        return;
    end
    fprintf('nplotKeogram: pooled %d IPP ROTI points.\n', numel(RO));

    tEdges = 0:1/6:(numel(dayRoots)*24);
    tC = (tEdges(1:end-1) + tEdges(2:end)) / 2;
    % zonal: lon bins
    lonEdges = floor(min(LO)):0.5:ceil(max(LO));
    lonC = (lonEdges(1:end-1) + lonEdges(2:end)) / 2;
    inB = AL >= latBand(1) & AL <= latBand(2);
    Z = binnedMax(TH(inB), LO(inB), RO(inB), tEdges, lonEdges);
    % meridional slices
    Ms = {};
    for q = 1:numel(lonSlices)
        inS = abs(LO - lonSlices(q)) <= 2.5;
        laE = floor(min(AL(inS))):0.5:ceil(max(AL(inS)));
        if numel(laE) < 3, laE = latBand(1):0.5:latBand(2); end
        laC = (laE(1:end-1) + laE(2:end)) / 2;
        Ms{q} = struct('M', binnedMax(TH(inS), AL(inS), RO(inS), tEdges, laE), ...
                       'ax', laC, 'lon', lonSlices(q)); %#ok<AGROW>
    end
    K = struct('tC', tC, 'lonC', lonC, 'Z', Z, 'slices', {Ms});

    out = fullfile(outDir, 'nplotKeogram.png');
    save(fullfile(outDir, 'keogram_data.mat'), '-struct', 'K');
    try
        fig = figure('Name', 'keogram', 'Position', [100 50 1300 350+260*(1+numel(lonSlices))], 'Visible', 'off');
        ax1 = subplot(1+numel(lonSlices), 1, 1);
        imagesc(tC, lonC, Z');
        set(ax1, 'YDir', 'normal');
        colorbar; caxis([0 1.5]); colormap(ax1, 'jet');
        ylabel('Longitude');
        title(sprintf('Zonal ROTI keogram, lat %.0f-%.0f (max/10 min / 0.5 deg)', latBand(1), latBand(2)));
        grid on;
        for q = 1:numel(lonSlices)
            ax = subplot(1+numel(lonSlices), 1, 1+q);
            imagesc(tC, Ms{q}.ax, Ms{q}.M');
            set(ax, 'YDir', 'normal');
            colorbar; caxis([0 1.5]); colormap(ax, 'jet');
            ylabel('Latitude');
            title(sprintf('Meridional slice near %gE', lonSlices(q)));
            grid on;
            if q == numel(lonSlices)
                xlabel('Hours past start (UT)');
            end
        end
        try
            saveas(fig, out);
        catch
            print(fig, out, '-dpng', '-r150', '-painters');
        end
        try, close(fig); catch, end
        fprintf('nplotKeogram saved: %s\n', out);
    catch ME
        fprintf('nplotKeogram figure failed (data saved): %s\n', ME.message);
        try, close all hidden; catch, end
    end
end

function Z = binnedMax(T, X, V, tEdges, xEdges)
% max V in each (time, x) bin; NaN where empty.
    Z = nan(numel(tEdges)-1, numel(xEdges)-1);
    [~, tb] = histc(T, tEdges);
    [~, xb] = histc(X, xEdges);
    ok = tb >= 1 & tb <= size(Z, 1) & xb >= 1 & xb <= size(Z, 2) & isfinite(V);
    tb = tb(ok); xb = xb(ok); V = V(ok);
    for k = 1:numel(V)
        if isnan(Z(tb(k), xb(k))) || V(k) > Z(tb(k), xb(k))
            Z(tb(k), xb(k)) = V(k);
        end
    end
end
