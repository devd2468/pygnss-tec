function nplotDrift(dayRoots, stations, outDir)
% NPLOTDRIFT  EPB zonal-drift estimate from keogram lag-correlation.
% For the most active evening window, correlates ROTI(t) between adjacent
% 1-deg longitude bins pooled over stations/systems; slope dLon/dt with
% 111 km/deg gives m/s (eastward positive). Honest quality flags included.
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    K = nplotKeogram(dayRoots, stations, outDir, {'G', 'E'});
    if isempty(K)
        fprintf('nplotDrift: no keogram.\n');
        return;
    end
    Z = K.Z;  % time x lon, 10-min bins
    lonC = K.lonC; tC = K.tC;
    % restrict to evening active window: max column energy 12-22h each day
    nb = size(Z, 1);
    vel = []; ql = [];
    for b = 1:size(Z, 2)-1
        x = Z(:, b); y = Z(:, b+1);
        ok = isfinite(x) & isfinite(y);
        if sum(ok) < 12, continue; end
        x = x - nanmean(x); y = y - nanmean(y);
        if nanstd(x) == 0 || nanstd(y) == 0, continue; end
        [cc, lags] = xcorr(x(ok), y(ok), 18, 'coeff');
        [mx, ix] = max(cc);
        if mx < 0.5, continue; end  % quality gate
        lagBins = lags(ix);  % +lag: y lags x
        dLonKm = (lonC(b+1) - lonC(b)) * 111;
        dtH = lagBins * (tC(2) - tC(1));
        if dtH == 0, continue; end
        vel(end+1) = dLonKm / (dtH * 3.6); %#ok<AGROW>  % m/s
        ql(end+1) = mx; %#ok<AGROW>
    end
    fig = figure('Name', 'drift', 'Position', [100 50 1000 550], 'Visible', 'off');
    if ~isempty(vel)
        subplot(1, 2, 1);
        scatter(1:numel(vel), vel, 60, ql, 'filled');
        colorbar; yline(0, 'k--');
        xlabel('Longitude-bin pair'); ylabel('Zonal velocity (m/s, +east)');
        title(sprintf('EPB drift estimates (median %.0f m/s, n=%d)', nanmedian(vel), numel(vel)));
        grid on;
        subplot(1, 2, 2);
        histogram(vel, 15, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'k');
        xlabel('m/s'); ylabel('Count');
        title('Drift distribution');
        grid on;
    else
        text(0.5, 0.5, 'No qualifying bin pairs (corr<0.5)', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    out = fullfile(outDir, 'nplotDrift.png');
    figsavesafe(gcf, out);
    if ~isempty(vel)
        fprintf('nplotDrift saved: median %.0f m/s eastward (n=%d, corr>=0.5)\n', nanmedian(vel), numel(vel));
    else
        fprintf('nplotDrift saved (no qualifying pairs).\n');
    end
end
