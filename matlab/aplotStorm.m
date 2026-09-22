function aplotStorm(dayRoots, stations, swMergedCsv, outDir)
% APLOTSTORM  Multi-day storm-context figure (GIX/ROTI/Dst style).
%   (a) VTEC median per station
%   (b) ROTI all-sat scatter + Dst (twin axis)
%   (c) Kp bars + IMF Bz
% Gray shading marks Dst < -50 nT intervals.
% Inputs:
%   dayRoots  - cell array of day folders (chronological)
%   stations  - cell array of 4-char codes (union across days)
%   swMergedCsv - space_weather_merged.csv path
%   outDir    - output folder
    if nargin < 4 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % ---- space weather ----
    [tSW, Kp]  = areadSW(swMergedCsv, 'Kp');
    [~, Dst]   = areadSW(swMergedCsv, 'Dst');
    [~, Bz]    = areadSW(swMergedCsv, 'IMF_Bz');
    [~, F107]  = areadSW(swMergedCsv, 'F10.7');
    if isempty(tSW)
        fprintf('aplotStorm: no space weather data in %s\n', swMergedCsv);
        return;
    end
    t0 = dateshift(min(tSW), 'start', 'day');
    swH = hours(tSW - t0);

    % ---- storm intervals (Dst < -50) ----
    shade = [];
    if ~isempty(Dst)
        bad = Dst < -50;
        bad(~isfinite(Dst)) = false;
        d = diff([false; bad; false]);
        st = find(d == 1); en = find(d == -1) - 1;
        for k = 1:numel(st)
            shade(end+1, :) = [swH(st(k)), swH(en(k))]; %#ok<AGROW>
        end
    end

    % ---- TEC/ROTI across days ----
    cmap = lines(max(numel(stations), 7));
    fig = figure('Name', 'Storm context', ...
        'Position', [100 50 1200 1000], 'Visible', 'off');

    % (a) VTEC medians
    ax1 = subplot(3, 1, 1); hold on; grid on;
    tEnd = 0;
    for d = 1:numel(dayRoots)
        for s = 1:numel(stations)
            S = aloadTEC(dayRoots{d}, stations{s});
            if isempty(S) || ~isfield(S, 'TEC'), continue; end
            v = nanmedian(S.TEC.vertical, 2);
            t = (d-1)*24 + (0:numel(v)-1)/3600;
            tEnd = max(tEnd, max(t));
            plot(t, v, '-', 'Color', cmap(s, :), 'LineWidth', 1.1);
        end
    end
    shadeX(ax1, shade);
    ylabel('VTEC (TECU)');
    title('Multi-day VTEC / ROTI with storm indices (GPS)');
    legend(stations, 'Location', 'best', 'FontSize', 7, 'NumColumns', 4);

    % (b) ROTI scatter + Dst
    ax2 = subplot(3, 1, 2); hold on; grid on;
    for d = 1:numel(dayRoots)
        for s = 1:numel(stations)
            S = aloadTEC(dayRoots{d}, stations{s});
            if isempty(S) || ~isfield(S, 'ROTI'), continue; end
            R = S.ROTI;
            keep = 1:30:size(R, 1);
            for p = 1:size(R, 2)
                v = R(keep, p);
                ok = ~isnan(v) & v > 0.02;
                if any(ok)
                    t = (d-1)*24 + (keep(ok)-1)/3600;
                    plot(t, v(ok), '.', 'Color', cmap(s, :), 'MarkerSize', 3);
                end
            end
        end
    end
    shadeX(ax2, shade);
    ylabel('ROTI (TECU/min)');
    if ~isempty(Dst)
        yyaxis right;
        plot(swH, Dst, 'k-', 'LineWidth', 1.4);
        ylabel('Dst (nT)');
        yyaxis left;
    end

    % (c) Kp bars + IMF Bz
    ax3 = subplot(3, 1, 3); hold on; grid on;
    if ~isempty(Kp)
        ok = isfinite(Kp);
        bar(swH(ok), Kp(ok), 1.5, 'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'none');
    end
    shadeX(ax3, shade);
    ylabel('Kp');
    if ~isempty(Bz)
        yyaxis right;
        plot(swH, Bz, 'r-', 'LineWidth', 1.1);
        ylabel('IMF Bz (nT)');
        yyaxis left;
    end
    xlabel('Hours past start (UT)');
    if ~isempty(F107)
        f0 = nanmedian(F107);
        title(ax1, sprintf('Multi-day VTEC / ROTI  (F10.7 median %.0f sfu)', f0));
    end

    linkaxes([ax1, ax2, ax3], 'x');
    xlim(ax1, [0, max(tEnd, max(swH))]);
    out = fullfile(outDir, 'aplotStorm_multiday.png');
    saveas(gcf, out);
    close;
    fprintf('aplotStorm saved: %s\n', out);
end

function shadeX(ax, shade)
    if isempty(shade), return; end
    yl = get(ax, 'YLim');
    if diff(yl) == 0, yl = [-1, 1]; end
    axes(ax); hold on;
    for k = 1:size(shade, 1)
        patch([shade(k,1) shade(k,2) shade(k,2) shade(k,1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.7 0.7 0.7], ...
              'EdgeColor', 'none', 'FaceAlpha', 0.35);
    end
    uistack(findobj(ax, 'Type', 'patch'), 'bottom');
    set(ax, 'YLim', yl);
end
