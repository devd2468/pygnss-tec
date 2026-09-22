function aplotKpBg(dayRoots, stations, swMergedCsv, outDir)
% APLOTKPBg  ROTI + VTEC over NOAA G-scale Kp background (storm-style).
% One stacked panel per station: ROTI dots (left axis) + VTEC line (right),
% background colored by 3-hourly Kp. Right colorbar shows G-scale.
    if nargin < 4 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    [tSW, Kp] = areadSW(swMergedCsv, 'Kp');
    if isempty(tSW)
        fprintf('aplotKpBg: no Kp data in %s\n', swMergedCsv);
        return;
    end
    t0 = dateshift(min(tSW), 'start', 'day');
    kpH = hours(tSW - t0);
    % Kp is 3-hourly: step edges
    nS = numel(stations);
    fig = figure('Name', 'Kp background', ...
        'Position', [100 50 1250, 250 + 260*nS], 'Visible', 'off');
    tEnd = max(kpH) + 3;
    for s = 1:nS
        ax = subplot(nS, 1, s); hold on;
        % background patches per Kp interval
        for k = 1:numel(kpH)
            if ~isfinite(Kp(k)), continue; end
            x0 = kpH(k); x1 = kpH(k) + 3;
            patch([x0 x1 x1 x0], [0 0 1 1], kpColor(Kp(k)), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.55);
        end
        % station data across days
        for d = 1:numel(dayRoots)
            S = aloadTEC(dayRoots{d}, stations{s});
            if isempty(S) || ~isfield(S, 'ROTI'), continue; end
            R = S.ROTI;
            keep = 1:20:size(R, 1);
            for p = 1:size(R, 2)
                v = R(keep, p);
                ok = ~isnan(v) & v > 0.02;
                if any(ok)
                    xx = (d-1)*24 + (keep(ok)-1)/3600;
                    plot(ax, xx, v(ok), 'b.', 'MarkerSize', 2);
                end
            end
            if isfield(S, 'TEC')
                vv = nanmedian(S.TEC.vertical, 2);
                tt = (d-1)*24 + (0:numel(vv)-1)/3600;
                yyaxis right;
                plot(ax, tt, vv, 'Color', [0 0.6 0.65], 'LineWidth', 1.4);
                ylabel(ax, 'VTEC');
                yyaxis left;
            end
        end
        % patches were drawn first: already behind data; no uistack
        % (uistack conflicts with yyaxis overlay axes)
        xlim([0, tEnd]);
        ylabel('ROTI');
        title(sprintf('%s', stations{s}));
        grid on;
        if s == nS
            xlabel('GPS time (hours past start)');
        end
    end
    % NOAA scale colorbar (discrete)
    colormap(ax, kpCmap());
    cb = colorbar('eastoutside');
    cb.Ticks = linspace(0, 1, 6);
    cb.TickLabels = {'Kp < 5', 'Kp = 5 (G1)', 'Kp = 6 (G2)', ...
        'Kp = 7 (G3)', 'Kp = 8 (G4)', 'Kp = 9 (G5)'};
    cb.Label.String = 'NOAA Scales Geomagnetic Storms';
    out = fullfile(outDir, 'aplotKpBg_stations.png');
    saveas(gcf, out);
    close;
    fprintf('aplotKpBg saved: %s\n', out);
end

function c = kpColor(kp)
% NOAA G-scale background colors: green -> yellow -> orange -> red -> darkred
    if kp < 5
        c = [0.55 0.9 0.35];
    elseif kp < 6
        c = [1.0 0.9 0.2];
    elseif kp < 7
        c = [1.0 0.65 0.1];
    elseif kp < 8
        c = [1.0 0.25 0.15];
    elseif kp < 9
        c = [0.75 0.1 0.1];
    else
        c = [0.45 0.0 0.0];
    end
end

function m = kpCmap()
    m = [0.55 0.9 0.35;
         1.0 0.9 0.2;
         1.0 0.65 0.1;
         1.0 0.25 0.15;
         0.75 0.1 0.1;
         0.45 0.0 0.0];
end
