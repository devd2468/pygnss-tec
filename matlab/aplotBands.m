function aplotBands(dayRoots, outDir)
% APLOTBANDS  Latitude-band gradient-style statistics (GIX/SIDX/ROTI layout).
% Stations are grouped by geographic latitude band; per 5-min bin:
%   col 1: band VTEC median (+/- IQR spread)
%   col 2: band ROTI median + ROTI 95th percentile
% Honest note: bands use station latitude (no per-IPP azimuth available).
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % ---- station latitudes (PPPindex next to this file; portable) ----
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
    sp = loadStationPositions(fullfile(here, 'PPPindex.txt'));
    if isempty(sp), sp = getDefaultStationData(); end
    stNames = {}; stLat = [];
    for k = 1:numel(sp)
        ll = ecef2lla(sp(k).xyz);
        stNames{end+1} = sp(k).name; %#ok<AGROW>
        stLat(end+1) = ll(1); %#ok<AGROW>
    end
    edges = [5, 15, 25, 35];
    nb = numel(edges) - 1;

    % ---- collect per-band series across days ----
    % bins: 5-min over full span
    nDays = numel(dayRoots);
    nBins = nDays * 288;
    Vmed = nan(nBins, nb); Vspr = nan(nBins, nb);
    Rmed = nan(nBins, nb); Rp95 = nan(nBins, nb);
    for d = 1:nDays
        files = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
        for f = 1:numel(files)
            try, S = load(fullfile(files(f).folder, files(f).name)); catch, continue; end
            if ~isfield(S, 'TEC'), continue; end
            st = S.station_name;
            ix = find(strcmp(stNames, st), 1);
            if isempty(ix), ix = find(strncmp(stNames, st, 4), 1); end
            if isempty(ix), continue; end
            b = find(stLat(ix) >= edges(1:end-1) & stLat(ix) < edges(2:end), 1);
            if isempty(b), continue; end
            V = S.TEC.vertical; R = S.ROTI;
            for k = 1:288
                i0 = (k-1)*300+1; i1 = min(k*300, size(V, 1));
                g = (d-1)*288 + k;
                vv = V(i0:i1, :); vv = vv(:); vv = vv(isfinite(vv));
                if ~isempty(vv)
                    Vmed(g, b) = nanmedian(vv);
                    Vspr(g, b) = iqr(vv);
                end
                rr = R(i0:min(i1, size(R, 1)), :); rr = rr(:); rr = rr(isfinite(rr));
                if ~isempty(rr)
                    Rmed(g, b) = nanmedian(rr);
                    Rp95(g, b) = prctile(rr, 95);
                end
            end
        end
    end
    tH = (0:nBins-1) * 5 / 60;

    fig = figure('Name', 'Band stats', ...
        'Position', [100 50 1300 350 + 260*nb], 'Visible', 'off');
    for b = 1:nb
        ax1 = subplot(nb, 2, 2*b-1); hold on; grid on;
        plot(tH, Vmed(:, b), 'k-', 'LineWidth', 1.2);
        plot(tH, Vmed(:, b) + Vspr(:, b)/2, 'r-', 'LineWidth', 0.8);
        plot(tH, max(Vmed(:, b) - Vspr(:, b)/2, 0), 'r-', 'LineWidth', 0.8);
        ylabel('VTEC (TECU)');
        title(sprintf('Band %g\\circ-%g\\circN  VTEC med/spread', edges(b), edges(b+1)));
        ax2 = subplot(nb, 2, 2*b); hold on; grid on;
        plot(tH, Rmed(:, b), 'g-', 'LineWidth', 1.2);
        plot(tH, Rp95(:, b), 'm-', 'LineWidth', 1.0);
        ylabel('ROTI (TECU/min)');
        title(sprintf('Band %g\\circ-%g\\circN  ROTI med/95', edges(b), edges(b+1)));
        if b == nb
            xlabel(ax1, 'Hours past start (UT)');
            xlabel(ax2, 'Hours past start (UT)');
        end
    end
    out = fullfile(outDir, 'aplotBands_multiday.png');
    saveas(gcf, out);
    close;
    fprintf('aplotBands saved: %s\n', out);
end
