function nplotROTImap10(dayRoots, stations, outDir, sysList, tSlots)
% NPLOTROTIMAP10  Abadi-style 10-min ROTI maps: 0.25-deg grid, max ROTI per
% cell, 5x5-cell smoothing, one panel per time slot.
% Defaults: sysList={'G'}, tSlots = 12:2:22 (hours past first-day 00 UT).
    if nargin < 4 || isempty(sysList), sysList = {'G'}; end
    if nargin < 5 || isempty(tSlots), tSlots = 12:2:22; end
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    AL = []; LO = []; TH = []; RO = [];
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
        fprintf('nplotROTImap10: no data.\n');
        return;
    end
    lonEdges = floor(min(LO)):0.25:ceil(max(LO));
    latEdges = floor(min(AL)):0.25:ceil(max(AL));
    lonC = (lonEdges(1:end-1) + lonEdges(2:end)) / 2;
    latC = (latEdges(1:end-1) + latEdges(2:end)) / 2;

    nP = numel(tSlots);
    nC = ceil(sqrt(nP)); nR = ceil(nP / nC);
    fig = figure('Name', 'rotimaps', 'Position', [100 50 380*nC+80, 320*nR+80], 'Visible', 'off');
    for k = 1:nP
        inT = TH >= tSlots(k) & TH < tSlots(k) + 1/6;
        M = binnedMax2(LO(inT), AL(inT), RO(inT), lonEdges, latEdges);
        % 5x5 smoothing (mean ignoring NaN)
        Ms = smooth5(M);
        ax = subplot(nR, nC, k); hold on;
        pcolor(lonC, latC, Ms');
        shading flat;
        colormap(ax, 'jet');
        caxis([0 1.5]);
        plot([96 108], [0 0], 'k-', 'LineWidth', 1);
        xlabel('Lon'); ylabel('Lat');
        title(sprintf('+%.1fh UT', tSlots(k)));
    end
    colorbar('eastoutside');
    sgtitle('10-min ROTI maps (max/0.25-deg cell, 5x5 smoothed)');
    out = fullfile(outDir, 'nplotROTImap10.png');
    figsavesafe(gcf, out);
    fprintf('nplotROTImap10 saved: %s\n', out);
end

function M = binnedMax2(X, Y, V, xE, yE)
    M = nan(numel(xE)-1, numel(yE)-1);
    if isempty(V), return; end
    [~, xb] = histc(X, xE);
    [~, yb] = histc(Y, yE);
    ok = xb >= 1 & xb <= size(M, 1) & yb >= 1 & yb <= size(M, 2) & isfinite(V);
    xb = xb(ok); yb = yb(ok); V = V(ok);
    for k = 1:numel(V)
        if isnan(M(xb(k), yb(k))) || V(k) > M(xb(k), yb(k))
            M(xb(k), yb(k)) = V(k);
        end
    end
end

function Ms = smooth5(M)
% 5x5 mean ignoring NaN (Abadi-style smoothing).
    Ms = M;
    [nx, ny] = size(M);
    Mp = nan(nx+4, ny+4);
    Mp(3:end-2, 3:end-2) = M;
    for i = 1:nx
        for j = 1:ny
            w = Mp(i:i+4, j:j+4);
            w = w(isfinite(w));
            if ~isempty(w)
                Ms(i, j) = mean(w);
            end
        end
    end
end
