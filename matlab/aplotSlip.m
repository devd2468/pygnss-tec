function aplotSlip(dayRoot, station4, outDir)
% APLOTSLIP  Cycle-slip proxy + visible-satellite count (PONC-style b/c).
%   Top:    per-5min cycle-slip counts (MW jump > 0.75 cyc) + NSat (twin axis)
%   Bottom: per-satellite ROTI dots
% Slip proxy uses Melbourne-Wuebbena jumps recomputed from saved obs.
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    S = aloadTEC(dayRoot, station4);
    if isempty(S) || ~isfield(S, 'obs')
        fprintf('aplotSlip: no obs for %s in %s\n', station4, dayRoot);
        return;
    end
    [GF, MW, prns, tag] = agf_mw(S.obs); %#ok<NASGU>
    if isempty(prns)
        fprintf('aplotSlip: no dual-frequency data for %s\n', station4);
        return;
    end
    obs = S.obs;
    % NSat per epoch from phase presence
    tL1 = pickFirst(obs.type, {'L1C', 'L1X'});
    L1c = obs.data(:, ismember(obs.type, tL1));
    sod = round(obs.epoch);
    epU = unique(sod);
    nsat = zeros(size(epU));
    for k = 1:numel(epU)
        nsat(k) = sum(~isnan(L1c(sod == epU(k))));
    end
    % slip counts per 5-min bin from MW jumps
    nB = 288;
    slips = zeros(nB, 1);
    for p = prns
        m = MW(:, p);
        iv = find(isfinite(m));
        if numel(iv) < 2, continue; end
        dm = abs(m(iv(2:end)) - m(iv(1:end-1)));
        ev = dm > 0.75;
        idx = iv(2:end);
        idx = idx(ev);
        b = min(nB, floor((idx-1)/300)+1);
        for k = 1:numel(b), slips(b(k)) = slips(b(k)) + 1; end
    end
    tB = (0:nB-1) * 5 / 60;

    fig = figure('Name', sprintf('Slips %s', station4), ...
        'Position', [100 50 1000 750], 'Visible', 'off');
    % top: slips + nsat
    ax1 = subplot(2, 1, 1); hold on; grid on;
    bar(tB, slips, 1, 'FaceColor', [0.5 0.3 0.7], 'EdgeColor', 'none');
    ylabel('Cycle Slips', 'Color', [0.5 0.3 0.7]);
    yyaxis right;
    stairs(epU/3600, nsat, 'Color', [0 0.7 0.7], 'LineWidth', 1.1);
    ylabel('NSat', 'Color', [0 0.7 0.7]);
    yyaxis left;
    title(sprintf('%s  cycle slips + visible sats  (%s)', station4, tag));
    xlim([0, 24]);
    % bottom: ROTI per sat
    ax2 = subplot(2, 1, 2); hold on; grid on;
    if isfield(S, 'ROTI')
        R = S.ROTI;
        cmap = lines(max(prns));
        keep = 1:10:size(R, 1);
        for p = prns
            v = R(keep, p);
            ok = ~isnan(v);
            if any(ok)
                plot((keep(ok)-1)/3600, v(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 4);
            end
        end
    end
    xlabel('Universal Time (hour)');
    ylabel('ROTI (TECU/min)');
    xlim([0, 24]);
    linkaxes([ax1, ax2], 'x');
    text(ax1, 0.02, 0.9, sprintf('Total slips: %d', sum(slips)), ...
        'Units', 'normalized', 'FontSize', 11, 'Color', [0.5 0.3 0.7]);

    out = fullfile(outDir, sprintf('aplotSlip_%s.png', station4));
    saveas(gcf, out);
    close;
    fprintf('aplotSlip saved: %s (total slips=%d)\n', out, sum(slips));
end

function t = pickFirst(types, cands)
    t = '';
    for k = 1:numel(cands)
        if any(strcmp(types, cands{k})), t = cands{k}; return; end
    end
end
