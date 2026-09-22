function aplotEvent(dayRoot, station4, tlimUTC, outDir)
% APLOTEVENT  HKWS/METU-style 3-panel event figure for one station-day.
%   Panel 1: per-satellite ROTI dots (TECU/min), gray 0-0.2 band
%   Panel 2: per-satellite ROTI^2 dots (TECU^2/min^2)
%   Panel 3: per-satellite DTEC lines (TECU), red box
% Inputs:
%   dayRoot  - e.g. 'D:\IGS DATA\May 2025\2025-05-01'
%   station4 - e.g. 'BHPL'
%   tlimUTC  - [h1 h2] window in hours (default full day dynamic range)
%   outDir   - output folder (default <dayRoot>/Results)
    if nargin < 3 || isempty(tlimUTC), tlimUTC = []; end
    if nargin < 4 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    S = aloadTEC(dayRoot, station4);
    if isempty(S) || ~isfield(S, 'TEC')
        fprintf('aplotEvent: no data for %s in %s\n', station4, dayRoot);
        return;
    end
    V = S.TEC.vertical;   % 86400 x 32 STEC-derived VTEC per PRN
    R = S.ROTI;           % 86400 x 32
    E = S.prm.elevation;  %#ok<NASGU>

    nPRN = size(V, 2);
    tH = (0:size(V, 1)-1) / 3600;
    cmap = lines(max(nPRN, 7));

    % per-satellite series, lightly decimated for dot plots
    keep = 1:30:size(V, 1);

    % auto window: centre on max ROTI if no window given
    if isempty(tlimUTC)
        rm = nanmedian(R, 2);
        [~, imx] = nanmax(rm);
        if isempty(imx) || isnan(imx), imx = 43200; end
        c = imx / 3600;
        tlimUTC = [max(0, c-1.5), min(24, c+1.5)];
        if diff(tlimUTC) < 1, tlimUTC = [0, 24]; end
    end
    inW = tH >= tlimUTC(1) & tH <= tlimUTC(2);

    fig = figure('Name', sprintf('Event %s', station4), ...
        'Position', [100 50 900 950], 'Visible', 'off');

    % ---- panel 1: ROTI per sat ----
    ax1 = subplot(3, 1, 1); hold on; grid on;
    for p = 1:nPRN
        v = R(keep, p);
        ok = ~isnan(v);
        if any(ok)
            plot(tH(keep(ok)), v(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 4);
        end
    end
    grayband(ax1, 0.2);
    ylabel('ROTI (TECU/min)');
    title(sprintf('%s  (GPS L1/L2, elev>30\\circ)', station4));
    xlim(tlimUTC);

    % ---- panel 2: ROTI^2 per sat ----
    ax2 = subplot(3, 1, 2); hold on; grid on;
    for p = 1:nPRN
        v = R(keep, p).^2;
        ok = ~isnan(v);
        if any(ok)
            plot(tH(keep(ok)), v(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 4);
        end
    end
    grayband(ax2, 0.04);
    ylabel('ROTI^2 (TECU^2/min^2)');

    % ---- panel 3: DTEC per sat ----
    ax3 = subplot(3, 1, 3); hold on; grid on;
    legN = {};
    for p = 1:nPRN
        vs = V(:, p);
        if sum(~isnan(vs)) < 60, continue; end
        d = adtec(vs, 3600);
        ok = ~isnan(d) & inW;
        if any(ok)
            plot(tH(ok), d(ok), '-', 'Color', cmap(p, :), 'LineWidth', 1.1);
            legN{end+1} = sprintf('G%02d', p); %#ok<AGROW>
        end
    end
    xlim(tlimUTC);
    xlabel('Universal Time (hour)');
    ylabel('DTEC (TECU)');
    % red box like the reference figures
    yl = get(ax3, 'YLim');
    rectangle('Position', [tlimUTC(1), yl(1), diff(tlimUTC), diff(yl)], ...
        'EdgeColor', 'r', 'LineWidth', 1.2);
    if ~isempty(legN)
        legend(legN, 'Location', 'best', 'FontSize', 7, 'NumColumns', 2);
    end
    % annotate strongest DTEC excursion as candidate TID/EPB
    try
        dd = adtec(nanmedian(V, 2), 3600);
        dd(~inW) = NaN;
        [mv, im] = nanmax(abs(dd));
        if ~isempty(mv) && ~isnan(mv) && mv > 0.25
            plot(tH(im), dd(im), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
            text(tH(im), dd(im), '  TIDs?', 'FontSize', 10, 'FontWeight', 'bold');
        end
    catch
    end

    linkaxes([ax1, ax2, ax3], 'x');
    xlim(ax1, tlimUTC);
    out = fullfile(outDir, sprintf('aplotEvent_%s_%02d-%02dUT.png', ...
        station4, round(tlimUTC(1)), round(tlimUTC(2))));
    saveas(gcf, out);
    close;
    fprintf('aplotEvent saved: %s\n', out);
end

function grayband(ax, y0)
% light-gray 0..y0 threshold band behind data
    xl = get(ax, 'XLim');
    if diff(xl) == 0, xl = [0, 24]; end
    patch(ax, [xl(1) xl(2) xl(2) xl(1)], [0 0 y0 y0], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    uistack(findobj(ax, 'Type', 'patch'), 'bottom');
end
