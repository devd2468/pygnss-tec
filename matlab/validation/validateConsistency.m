function validateConsistency(dayRoots, outDir)
% VALIDATECONSISTENCY  Local-only validation suite (no internet required).
% Four GPS-Solutions-grade internal checks using only local products:
%   1. Receiver-DCB day-to-day stability per station (same hardware ->
%      stable estimates; spread quantifies estimator uncertainty).
%   2. Inter-station VTEC agreement: pairwise correlation of 5-min median
%      VTEC series per day (nearby stations must track, r > ~0.95).
%   3. Before/after OSB: VTEC level shift from dcb_before.csv (pre-OSB run).
%   4. Constellation agreement: Galileo/BeiDou/QZSS VTEC vs GPS VTEC per
%      station-day (Pearson r, median bias, RMSE on overlapping 5-min bins).
% Writes validation_stats.csv + validation_consistency.png.
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % ---- gather per station-day (GPS base files only; SYS files handled in §4) ----
    recs = struct('day', {}, 'station', {}, 'vmean', {}, 'rcv', {}, 'vmed', {});
    for d = 1:numel(dayRoots)
        dd = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
        for k = 1:numel(dd)
            if ~isempty(regexp(dd(k).name, '^MultiGNSS_[A-Z0-9]{4}_\d{4}_', 'once'))
                % GPS base file only (skip MultiGNSS_<ST>_<SYS>_*.mat here)
            else
                continue;
            end
            try
                S = load(fullfile(dd(k).folder, dd(k).name));
                V = S.TEC.vertical;
                vm = fiveMinMed(V);
                recs(end+1) = struct('day', d, 'station', S.station_name, ...
                    'vmean', nanmean(V(:)), 'rcv', S.DCB.rcv, 'vmed', vm); %#ok<AGROW>
            catch
            end
        end
    end
    if isempty(recs)
        fprintf('validateConsistency: no data.\n');
        return;
    end

    sts = unique({recs.station});
    % ---- 1. rcv stability ----
    fprintf('\nReceiver-DCB day-to-day stability (TECU):\n');
    stab = nan(numel(sts), 1);
    for s = 1:numel(sts)
        ix = find(strcmp({recs.station}, sts{s}));
        vals = [recs(ix).rcv];
        stab(s) = std(vals);
        fprintf('  %-5s n=%d mean=%7.2f std=%5.2f range=[%.1f, %.1f]\n', ...
            sts{s}, numel(vals), mean(vals), stab(s), min(vals), max(vals));
    end

    % ---- 2. inter-station agreement (per day, pairwise r of 5-min medians) ----
    fprintf('\nInter-station VTEC agreement (Pearson r of 5-min medians):\n');
    agree = [];
    for d = 1:numel(dayRoots)
        ix = find([recs.day] == d);
        for a = 1:numel(ix)
            for b = a+1:numel(ix)
                x = recs(ix(a)).vmed; y = recs(ix(b)).vmed;
                ok = isfinite(x) & isfinite(y);
                if sum(ok) > 60
                    C = corrcoef(x(ok), y(ok));
                    agree(end+1, :) = [d, C(1, 2)]; %#ok<AGROW>
                end
            end
        end
    end
    if ~isempty(agree)
        for d = 1:numel(dayRoots)
            r = agree(agree(:, 1) == d, 2);
            if ~isempty(r)
                fprintf('  day%d: pairs=%d median r=%.4f min r=%.4f\n', ...
                    d, numel(r), nanmedian(r), nanmin(r));
            end
        end
    end

    % ---- 3. before/after OSB VTEC shift ----
    bvFile = fullfile(dayRoots{1}, '..', 'ArticleFigs', 'dcb_before.csv');
    dV = [];
    if exist(bvFile, 'file')
        try
            T = readtable(bvFile, 'TextType', 'string');
            for k = 1:numel(recs)
                m = find(T.day == recs(k).day & strcmp(string(T.station), recs(k).station));
                if ~isempty(m)
                    dV(end+1) = recs(k).vmean - T.vtec_mean(m(1)); %#ok<AGROW>
                end
            end
            fprintf('\nOSB effect on VTEC mean: n=%d median dV=%+.2f TECU\n', ...
                numel(dV), nanmedian(dV));
        catch ME
            fprintf('\nBefore/after compare skipped: %s\n', ME.message);
        end
    end

    % ---- 4. constellation agreement: E/C/J VTEC vs GPS VTEC ----
    % Per station-day, 5-min medians on overlapping bins: Pearson r,
    % median bias (SYS-GPS) and RMSE. Independent orbit/frequency chains
    % must track GPS (r > ~0.9); bias absorbs DCB-leveling differences.
    fprintf('\nConstellation agreement vs GPS (5-min VTEC medians):\n');
    sysList = {'E', 'C', 'J'};
    sysAgree = struct();
    for t = 1:numel(sysList)
        sysAgree.(sysList{t}).r = [];
        sysAgree.(sysList{t}).bias = [];
        sysAgree.(sysList{t}).rmse = [];
        sysAgree.(sysList{t}).n = 0;
    end
    for d = 1:numel(dayRoots)
        gd = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
        % GPS base files in this day
        for k = 1:numel(gd)
            m = regexp(gd(k).name, '^MultiGNSS_([A-Z0-9]{4})_(\d{4}_\d{2}_\d{2})\.mat$', 'tokens', 'once');
            if isempty(m), continue; end
            st = m{1}; ds = m{2};
            try
                G = load(fullfile(gd(k).folder, gd(k).name));
                gv = fiveMinMed(G.TEC.vertical);
            catch
                continue;
            end
            for t = 1:numel(sysList)
                sf = dir(fullfile(dayRoots{d}, 'Results', ...
                    sprintf('MultiGNSS_%s_%s_%s.mat', st, sysList{t}, ds)));
                if isempty(sf), continue; end
                try
                    S2 = load(fullfile(sf(1).folder, sf(1).name));
                    sv = fiveMinMed(S2.TEC.vertical);
                    n = min(numel(gv), numel(sv));
                    ok = isfinite(gv(1:n)) & isfinite(sv(1:n));
                    if sum(ok) < 60, continue; end
                    x = gv(1:n); x = x(ok); y = sv(1:n); y = y(ok);
                    C = corrcoef(x, y);
                    sysAgree.(sysList{t}).r(end+1) = C(1, 2); %#ok<AGROW>
                    sysAgree.(sysList{t}).bias(end+1) = nanmedian(y - x); %#ok<AGROW>
                    sysAgree.(sysList{t}).rmse(end+1) = sqrt(nanmean((y - x).^2)); %#ok<AGROW>
                    sysAgree.(sysList{t}).n = sysAgree.(sysList{t}).n + 1;
                catch
                end
            end
        end
    end
    for t = 1:numel(sysList)
        s = sysList{t};
        if sysAgree.(s).n > 0
            fprintf('  %s vs GPS: pairs=%d median r=%.4f | bias=%+.2f TECU | RMSE=%.2f TECU\n', ...
                s, sysAgree.(s).n, nanmedian(sysAgree.(s).r), ...
                nanmedian(sysAgree.(s).bias), nanmedian(sysAgree.(s).rmse));
        else
            fprintf('  %s vs GPS: no overlapping pairs.\n', s);
        end
    end

    % ---- save stats ----
    fid = fopen(fullfile(outDir, 'validation_stats.csv'), 'w');
    fprintf(fid, 'metric,station_or_day,value\n');
    for s = 1:numel(sts)
        fprintf(fid, 'rcv_std,%s,%.3f\n', sts{s}, stab(s));
    end
    for d = 1:numel(dayRoots)
        r = agree(agree(:, 1) == d, 2);
        if ~isempty(r)
            fprintf(fid, 'agree_median_r,day%d,%.4f\n', d, nanmedian(r));
            fprintf(fid, 'agree_min_r,day%d,%.4f\n', d, nanmin(r));
        end
    end
    if ~isempty(dV)
        fprintf(fid, 'osb_dV_median,all,%.3f\n', nanmedian(dV));
    end
    for t = 1:numel(sysList)
        s = sysList{t};
        if sysAgree.(s).n > 0
            fprintf(fid, 'sys_median_r,%s,%.4f\n', s, nanmedian(sysAgree.(s).r));
            fprintf(fid, 'sys_median_bias,%s,%.3f\n', s, nanmedian(sysAgree.(s).bias));
            fprintf(fid, 'sys_median_rmse,%s,%.3f\n', s, nanmedian(sysAgree.(s).rmse));
            fprintf(fid, 'sys_pairs,%s,%d\n', s, sysAgree.(s).n);
        end
    end
    fclose(fid);

    % ---- figure ----
    fig = figure('Name', 'validation', 'Position', [100 50 1200 850], 'Visible', 'off');
    subplot(2, 2, 1); hold on; grid on;
    for s = 1:numel(sts)
        ix = find(strcmp({recs.station}, sts{s}));
        plot([recs(ix).day], [recs(ix).rcv], 'o-', 'LineWidth', 1.2);
    end
    xlabel('Day (May 2025)'); ylabel('Receiver DCB (TECU)');
    title('Receiver-DCB day-to-day stability');
    legend(sts, 'Location', 'best', 'FontSize', 7);
    subplot(2, 2, 2); hold on; grid on;
    if ~isempty(agree)
        for d = 1:numel(dayRoots)
            r = agree(agree(:, 1) == d, 2);
            if ~isempty(r)
                plot(d*ones(size(r)), r, 'ko', 'MarkerSize', 5);
            end
        end
    end
    yline(0.95, 'r--', 'r=0.95');
    xlabel('Day'); ylabel('Pairwise r (VTEC medians)');
    title('Inter-station VTEC agreement');
    xlim([0.5, numel(dayRoots)+0.5]);
    subplot(2, 2, 3); hold on; grid on;
    if ~isempty(dV)
        histogram(dV, 15, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'k');
        xline(nanmedian(dV), 'r-', 'LineWidth', 2);
    end
    xlabel('\DeltaVTEC mean, after-before OSB (TECU)');
    ylabel('Station-days');
    title('OSB satellite-DCB effect on absolute VTEC');
    subplot(2, 2, 4); hold on; grid on;
    bn = {'Galileo (E)', 'BeiDou (C)', 'QZSS (J)'};
    for t = 1:numel(sysList)
        s = sysList{t};
        if sysAgree.(s).n > 0
            r = sysAgree.(s).r;
            bar(t, nanmedian(r), 'FaceColor', [0.2 0.6 0.4]);
            plot(t*ones(size(r)), r, 'ko', 'MarkerSize', 4);
            text(t, nanmedian(r)+0.01, ...
                sprintf('bias %+.1f\nRMSE %.1f', ...
                nanmedian(sysAgree.(s).bias), nanmedian(sysAgree.(s).rmse)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
    set(gca, 'XTick', 1:numel(sysList), 'XTickLabel', bn);
    yline(0.9, 'r--', 'r=0.9');
    ylabel('Median r vs GPS');
    title('Constellation VTEC agreement');
    ylim([0, 1.05]);
    out = fullfile(outDir, 'validation_consistency.png');
    saveas(gcf, out);
    close;
    fprintf('validation saved: %s\n', out);
end

function vm = fiveMinMed(V)
% 5-minute medians of an 86400xN VTEC matrix (column vector output).
    nb = floor(size(V, 1) / 300);
    vm = nan(nb, 1);
    for b = 1:nb
        blk = V((b-1)*300+1:min(b*300, size(V, 1)), :);
        vm(b) = nanmedian(blk(:));
    end
end
