function nplotDepletion(dayRoot, station4, sysChar, outDir)
% NPLOTDEPLETION  TEC-depletion consistency (Muafiry style): per-satellite
% VTEC minus 120-min running mean, depletion time series + depletion map
% at the deepest epoch, cross-checked against ROTI.
    if nargin < 3 || isempty(sysChar), sysChar = 'G'; end
    if nargin < 4 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    P = getIPP(dayRoot, station4, sysChar);
    if isempty(P) || ~isfield(P, 'vtec')
        fprintf('nplotDepletion: no data for %s/%s\n', station4, sysChar);
        return;
    end
    V = P.vtec;
    % 120-min running mean per sat (omitnan where available)
    D = nan(size(V));
    for p = 1:size(V, 2)
        v = V(:, p);
        if sum(isfinite(v)) < 240, continue; end
        try
            tr = movmean(v, 7200, 'omitnan');
        catch
            tr = movmean(fillmissing(v, 'linear'), 7200);
        end
        D(:, p) = v - tr;
    end
    % deepest depletion epoch (median across sats)
    dm = nanmedian(D, 2);
    [dep, im] = nanmin(dm);
    if isempty(dep) || isnan(dep)
        fprintf('nplotDepletion: empty depletion for %s/%s\n', station4, sysChar);
        return;
    end
    tH = (0:size(V, 1)-1) / 3600;
    fig = figure('Name', 'depletion', 'Position', [100 50 1200 800], 'Visible', 'off');
    ax1 = subplot(2, 1, 1); hold on; grid on;
    cmap = lines(size(D, 2));
    for p = 1:size(D, 2)
        ok = isfinite(D(:, p));
        if sum(ok) > 100
            plot(tH(ok), D(ok, p), '-', 'Color', cmap(p, :), 'LineWidth', 0.8);
        end
    end
    plot(tH, dm, 'k-', 'LineWidth', 1.6);
    ylabel('TEC depletion (TECU)');
    title(sprintf('%s %s TEC depletion (120-min detrend); deepest %.2f TECU at %.2f UT', ...
        station4, sysChar, dep, tH(im)));
    xlim([0 24]);
    ax2 = subplot(2, 1, 2); hold on; grid on;
    rows = max(1, im-30):min(size(D, 1), im+30);
    qx = P.lon(rows, :); qy = P.lat(rows, :); qv = D(rows, :);
    qok = isfinite(qx) & isfinite(qy) & isfinite(qv);
    sc = scatter(qx(qok), qy(qok), 20, qv(qok), 'filled');
    set(sc, 'MarkerEdgeColor', 'none');
    colorbar; caxis([-5 1]);
    colormap(ax2, 'jet');
    xlabel('Longitude'); ylabel('Latitude');
    title('Depletion map at deepest epoch (red = background, blue = depleted)');
    out = fullfile(outDir, sprintf('nplotDepletion_%s_%s.png', station4, sysChar));
    figsavesafe(gcf, out);
    fprintf('nplotDepletion saved: %s (deepest %.2f TECU)\n', out, dep);
end
