function aplotGFstats(dayRoot, station4, outDir)
% APLOTGFSTATS  Geometry-free residual statistics (article style):
%   Fig A: histogram of epoch-differenced GF phase (m) with Mean/STD.
%   Fig B: scatter of dPhi_GF (m) vs ROTI with red envelope bounds.
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    S = aloadTEC(dayRoot, station4);
    if isempty(S) || ~isfield(S, 'obs') || ~isfield(S, 'ROTI')
        fprintf('aplotGFstats: no data for %s in %s\n', station4, dayRoot);
        return;
    end
    [GF, ~, prns, tag] = agf_mw(S.obs);
    if isempty(prns)
        fprintf('aplotGFstats: no dual-frequency data for %s\n', station4);
        return;
    end
    R = S.ROTI;

    % epoch-differenced GF pooled over sats (consecutive VALID samples).
    % ROTI lives on an interleaved row grid -> nearest valid ROTI <= 120 s.
    dAll = []; rAll = [];
    for p = prns
        g = GF(:, p);
        iv = find(isfinite(g));
        if numel(iv) < 2, continue; end
        dg = g(iv(2:end)) - g(iv(1:end-1));
        iq = iv(2:end);
        rv = R(:, p);
        vidx = find(isfinite(rv));
        if isempty(vidx), continue; end
        rr = interp1(double(vidx), rv(vidx), double(iq), 'nearest', NaN);
        near = interp1(double(vidx), double(vidx), double(iq), 'nearest', NaN);
        ok = abs(dg) < 1.5 & isfinite(dg) & isfinite(rr) & abs(iq - near) <= 120;
        dAll = [dAll; dg(ok)]; %#ok<AGROW>
        rAll = [rAll; rr(ok)]; %#ok<AGROW>
    end
    ok = isfinite(dAll) & isfinite(rAll);
    dAll = dAll(ok); rAll = rAll(ok);
    if isempty(dAll)
        fprintf('aplotGFstats: empty residuals for %s\n', station4);
        return;
    end

    % ---- Fig A: histogram ----
    fig = figure('Name', 'GF hist', 'Position', [100 100 700 600], 'Visible', 'off');
    edges = -0.5:0.025:0.5;
    h = histogram(dAll, edges, 'FaceColor', [0 0 0.8], 'EdgeColor', 'k');
    % percentage axis
    yt = get(gca, 'YTick');
    set(gca, 'YTickLabel', arrayfun(@(y) sprintf('%.0f', y/numel(dAll)*100), yt, 'UniformOutput', false));
    mu = mean(dAll); sg = std(dAll);
    text(0.05, 0.9, sprintf('Mean = %.3f\nSTD = %.3f', mu, sg), ...
        'Units', 'normalized', 'FontSize', 16);
    xlabel('\Delta\Phi_{GF} (m)');
    ylabel('Percentage (%)');
    grid on;
    out = fullfile(outDir, sprintf('aplotGFhist_%s.png', station4));
    saveas(gcf, out);
    close;

    % ---- Fig B: scatter vs ROTI with red envelope ----
    fig = figure('Name', 'GF vs ROTI', 'Position', [100 100 750 650], 'Visible', 'off');
    hold on; grid on;
    % decimate dense cloud
    idx = randperm(numel(dAll));
    idx = idx(1:min(60000, numel(idx)));
    plot(rAll(idx), dAll(idx), 'b.', 'MarkerSize', 2);
    % red envelope: +/-0.1 flat to 0.5, ramp to +/-0.6 at 3, flat after
    xe = [0, 0.5, 3, 6];
    ye = [0.1, 0.1, 0.6, 0.6];
    plot(xe, ye, 'r-', 'LineWidth', 2);
    plot(xe, -ye, 'r-', 'LineWidth', 2);
    xlabel('ROTI (TECU/min)');
    ylabel('\Delta\Phi_{GF} (m)');
    title(sprintf('%s  GF residuals vs ROTI  (%s)', station4, tag));
    xlim([0, 6]); ylim([-1.5, 1.5]);
    out = fullfile(outDir, sprintf('aplotGFscatter_%s.png', station4));
    saveas(gcf, out);
    close;
    fprintf('aplotGFstats saved for %s (n=%d, mean=%.4f, std=%.4f)\n', ...
        station4, numel(dAll), mu, sg);
end
