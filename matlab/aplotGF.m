function aplotGF(dayRoot, station4, outDir)
% APLOTGF  SCOR-style 3-panel per-satellite figure:
%   (1) Melbourne-Wuebbena delta-N (cycles), ylim [-2,2]
%   (2) geometry-free delta-Phi (m), +/-0.05 m lines
%   (3) ROTI per satellite (TECU/min), PRN legend
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    S = aloadTEC(dayRoot, station4);
    if isempty(S) || ~isfield(S, 'obs')
        fprintf('aplotGF: no obs for %s in %s\n', station4, dayRoot);
        return;
    end
    [GF, MW, prns, tag] = agf_mw(S.obs);
    if isempty(prns)
        fprintf('aplotGF: no dual-frequency data for %s\n', station4);
        return;
    end
    R = [];
    if isfield(S, 'ROTI'), R = S.ROTI; end
    tH = (0:86399) / 3600;
    cmap = lines(max(prns));

    fig = figure('Name', sprintf('GF %s', station4), ...
        'Position', [100 50 1100 950], 'Visible', 'off');

    % (1) dN_MW (diffs of consecutive VALID samples: data are 30 s sampled)
    ax1 = subplot(3, 1, 1); hold on; grid on;
    for p = prns
        m = MW(:, p);
        iv = find(isfinite(m));
        if numel(iv) > 1
            dm = m(iv(2:end)) - m(iv(1:end-1));
            tv = tH(iv(2:end));
            big = abs(dm) > 2; dm(big) = NaN; tv(big) = NaN;  % bound slip spikes
            ok = ~isnan(dm);
            if any(ok)
                plot(tv(ok), dm(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 3);
            end
        end
    end
    ylabel('\DeltaN_{MW} (cycle)');
    title(sprintf('%s  MW / GF / ROTI  (%s)', station4, tag));
    ylim([-2, 2]); xlim([0, 24]);

    % (2) dPhi_GF (diffs of consecutive VALID samples)
    ax2 = subplot(3, 1, 2); hold on; grid on;
    for p = prns
        g = GF(:, p);
        iv = find(isfinite(g));
        if numel(iv) > 1
            dg = g(iv(2:end)) - g(iv(1:end-1));
            tv = tH(iv(2:end));
            big = abs(dg) > 0.6; dg(big) = NaN; tv(big) = NaN;
            ok = ~isnan(dg);
            if any(ok)
                plot(tv(ok), dg(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 3);
            end
        end
    end
    yline(0.05, 'k-'); yline(-0.05, 'k-');
    text(12, 0.3, '0.05 m', 'FontSize', 9);
    text(12, -0.3, '-0.05 m', 'FontSize', 9);
    ylabel('\Delta\\Phi_{GF} (m)');
    ylim([-0.6, 0.6]); xlim([0, 24]);

    % (3) ROTI per sat
    ax3 = subplot(3, 1, 3); hold on; grid on;
    legN = {};
    if ~isempty(R)
        keep = 1:10:size(R, 1);
        for p = prns
            v = R(keep, p);
            ok = ~isnan(v);
            if any(ok)
                plot(tH(keep(ok)), v(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 4);
                legN{end+1} = sprintf('G%02d', p); %#ok<AGROW>
            end
        end
    end
    ylabel('ROTI (TECU/min)');
    xlabel('GPS time (hour)');
    xlim([0, 24]);
    if ~isempty(legN)
        legend(legN, 'Location', 'eastoutside', 'FontSize', 6, 'NumColumns', 2);
    end
    linkaxes([ax1, ax2, ax3], 'x');

    out = fullfile(outDir, sprintf('aplotGF_%s.png', station4));
    saveas(gcf, out);
    close;
    fprintf('aplotGF saved: %s\n', out);
end
