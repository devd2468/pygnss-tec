function R = dcbWelchTests(dayRoots, outDir)
% DCBWELCHTESTS  Formal receiver-DCB / VTEC hypothesis tests for the paper.
% No toolbox dependency: t p-values via the incomplete beta function,
% plus fixed-seed permutation p-values (assumption-free robustness).
%
% Test A (paired): VTEC mean before vs after OSB satellite DCB (n=34 pairs).
%   H0: OSB correction does not shift absolute VTEC.
% Test B (Welch two-sample): receiver DCB on the active EPB night (04-May)
%   vs quiet days (H0: equal means -> Ma-Maruyama robust across conditions).
% Test C (paired per system): receiver DCB GPS vs Galileo/BeiDou/QZSS on
%   the same station-days, i.e. same hardware (H0: zero mean difference).
% Outputs: dcb_welch.txt (full report) + dcb_welch.png (forest plot w/ CI).
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    R = struct('test', {}, 'n', {}, 'est', {}, 'ci', {}, ...
               't', {}, 'df', {}, 'p', {}, 'pPerm', {}, 'd', {});
    % ---------------- Test A: paired VTEC before/after OSB ----------------
    bvFile = fullfile(dayRoots{1}, '..', 'ArticleFigs', 'dcb_before.csv');
    if exist(bvFile, 'file')
        try
            T = readtable(bvFile, 'TextType', 'string');
            after = []; before = [];
            for d = 1:numel(dayRoots)
                dd = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
                for k = 1:numel(dd)
                    if isempty(regexp(dd(k).name, '^MultiGNSS_[A-Z0-9]{4}_\d{4}_', 'once'))
                        continue;  % GPS base only
                    end
                    S = load(fullfile(dd(k).folder, dd(k).name));
                    m = find(T.day == d & strcmp(string(T.station), S.station_name));
                    if ~isempty(m)
                        before(end+1) = T.vtec_mean(m(1)); %#ok<AGROW>
                        after(end+1) = nanmean(S.TEC.vertical(:)); %#ok<AGROW>
                    end
                end
            end
            ddif = after - before;
            ddif = ddif(isfinite(ddif));
            if numel(ddif) >= 5
                [t, df, p, ci, dz, pp] = pairedTest(ddif);
                R(end+1) = pack('A: VTEC after-before OSB (paired)', ...
                    numel(ddif), mean(ddif), ci, t, df, p, pp, dz); %#ok<AGROW>
            end
        catch ME
            fprintf('Test A skipped: %s\n', ME.message);
        end
    else
        fprintf('Test A skipped: dcb_before.csv not found.\n');
    end

    % ---------------- Test B: Welch active-night vs quiet-days rcv DCB ----
    try
        act = []; qui = [];
        for d = 1:numel(dayRoots)
            dd = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
            for k = 1:numel(dd)
                if isempty(regexp(dd(k).name, '^MultiGNSS_[A-Z0-9]{4}_\d{4}_', 'once'))
                    continue;
                end
                S = load(fullfile(dd(k).folder, dd(k).name));
                v = S.DCB.rcv;
                if ~isscalar(v) || ~isfinite(v), continue; end
                if d == 4
                    act(end+1) = v; %#ok<AGROW>
                else
                    qui(end+1) = v; %#ok<AGROW>
                end
            end
        end
        if numel(act) >= 3 && numel(qui) >= 3
            [t, df, p, ci, d, pp] = welchTest(act, qui);
            R(end+1) = pack('B: rcvDCB active-night vs quiet (Welch)', ...
                numel(act)+numel(qui), mean(act)-mean(qui), ci, t, df, p, pp, d); %#ok<AGROW>
        end
    catch ME
        fprintf('Test B skipped: %s\n', ME.message);
    end

    % ---------------- Test C: paired GPS vs E/C receiver DCB ------------
    % QZSS handled separately below (Test C-J): its uniformly high-elevation
    % IGSO geometry flattens the Ma-Maruyama cost, so single-day searches
    % peg at the bounds; 5-day pooling recovers interior minima where the
    % geometry admits it (see cost-curve diagnostic).
    for s = {'E', 'C'}
        sys = s{1};
        dg = []; ds = [];
        for d = 1:numel(dayRoots)
            gd = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
            for k = 1:numel(gd)
                m = regexp(gd(k).name, '^MultiGNSS_([A-Z0-9]{4})_(\d{4}_\d{2}_\d{2})\.mat$', 'tokens', 'once');
                if isempty(m), continue; end
                st = m{1}; ds2 = m{2};
                sf = dir(fullfile(dayRoots{d}, 'Results', ...
                    sprintf('MultiGNSS_%s_%s_%s.mat', st, sys, ds2)));
                if isempty(sf), continue; end
                try
                    G = load(fullfile(gd(k).folder, gd(k).name));
                    S2 = load(fullfile(sf(1).folder, sf(1).name));
                    gv = G.DCB.rcv; sv = S2.DCB.rcv;
                    if isscalar(gv) && isscalar(sv) && isfinite(gv) && isfinite(sv)
                        dg(end+1) = gv; %#ok<AGROW>
                        ds(end+1) = sv; %#ok<AGROW>
                    end
                catch
                end
            end
        end
        if numel(dg) >= 5
            ddif = ds - dg;
            [t, df, p, ci, dz, pp] = pairedTest(ddif);
            R(end+1) = pack(sprintf('C: rcvDCB %s-GPS paired (same hardware)', sys), ...
                numel(ddif), mean(ddif), ci, t, df, p, pp, dz); %#ok<AGROW>
        else
            fprintf('Test C-%s skipped: only %d pairs.\n', sys, numel(dg));
        end
    end

    % -------- Test C-J: pooled-QZSS vs GPS receiver DCB (station-level) ---
    % Single-day J estimates peg at the search bounds (sparse 3-4 sat arcs),
    % so J is estimated once per station from all 5 days pooled (same
    % hardware assumption as GPS day-means). J TEC/ROTI products unchanged.
    try
        jg = []; js = [];
        % station list from GPS base files
        stAll = {};
        for d = 1:numel(dayRoots)
            gd = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
            for k = 1:numel(gd)
                m = regexp(gd(k).name, '^MultiGNSS_([A-Z0-9]{4})_\d{4}_', 'tokens', 'once');
                if ~isempty(m) && ~any(strcmp(stAll, m{1}))
                    stAll{end+1} = m{1}; %#ok<AGROW>
                end
            end
        end
        jst = {};
        for s = 1:numel(stAll)
            [rj, ~, ~, peg] = poolReceiverDCB(dayRoots, stAll{s}, 'J');
            if ~isfinite(rj) || peg, continue; end
            % GPS station mean across days
            gv = [];
            for d = 1:numel(dayRoots)
                gd = dir(fullfile(dayRoots{d}, 'Results', ...
                    sprintf('MultiGNSS_%s_2*.mat', stAll{s})));
                for k = 1:numel(gd)
                    if ~isempty(regexp(gd(k).name, '^MultiGNSS_[A-Z0-9]{4}_[ERCJ]_', 'once'))
                        continue;
                    end
                    try
                        G = load(fullfile(gd(k).folder, gd(k).name));
                        if isscalar(G.DCB.rcv) && isfinite(G.DCB.rcv)
                            gv(end+1) = G.DCB.rcv; %#ok<AGROW>
                        end
                    catch
                    end
                end
            end
            if ~isempty(gv)
                jg(end+1) = mean(gv); %#ok<AGROW>
                js(end+1) = rj; %#ok<AGROW>
                jst{end+1} = stAll{s}; %#ok<AGROW>
            end
        end
        % Descriptive pairs regardless of test eligibility (no p-value when n<4:
        % permutation floor makes formal testing vacuous; values stay honest).
        for q = 1:numel(jg)
            fprintf('  C-J %s: GPS mean %7.2f -> pooled J %7.2f (d=%+7.2f)\n', ...
                jst{q}, jg(q), js(q), js(q)-jg(q));
        end
        if numel(jg) >= 4
            ddif = js - jg;
            [t, df, p, ci, dz, pp] = pairedTest(ddif);
            R(end+1) = pack('C: rcvDCB J-GPS paired pooled (same hardware)', ...
                numel(ddif), mean(ddif), ci, t, df, p, pp, dz); %#ok<AGROW>
        else
            fprintf('Test C-J: n=%d below formal-test floor (permutation p floor 0.25); reported descriptively only.\n', numel(jg));
        end
    catch ME
        fprintf('Test C-J skipped: %s\n', ME.message);
    end

    % -------- diagnostic: cost curves, stabilized vs structural cases ----
    % Left (LCK4): single-day J unstable -> pooled interior minimum (fix works).
    % Right (BHPL): single-day AND pooled J pegged at bound (uniform slant
    % factors flatten the cost: structural identifiability limit).
    try
        fig2 = figure('Name', 'costcurves', 'Position', [100 50 1400 550], 'Visible', 'off');
        hn = @(c) (c - nanmin(c)) / max(nanmax(c) - nanmin(c), eps);
        demo = {'LCK4', 'BHPL'};
        subt = {'LCK4: pooling stabilizes the estimate', ...
                'BHPL: pegged single-day AND pooled (structural limit)'};
        for pp = 1:2
            subplot(1, 2, pp); hold on; grid on;
            [~, cG, gG] = poolReceiverDCB(dayRoots(1), demo{pp}, 'G');
            [~, cJ1, gJ1] = poolReceiverDCB(dayRoots(1), demo{pp}, 'J');
            [~, cJP, gJP] = poolReceiverDCB(dayRoots, demo{pp}, 'J');
            lg = {};
            if any(isfinite(cG))
                plot(gG, hn(cG), 'k-', 'LineWidth', 1.5);
                lg{end+1} = 'GPS single day'; %#ok<AGROW>
            end
            if any(isfinite(cJ1))
                plot(gJ1, hn(cJ1), 'r--', 'LineWidth', 1.2);
                lg{end+1} = 'QZSS single day'; %#ok<AGROW>
            end
            if any(isfinite(cJP))
                plot(gJP, hn(cJP), 'b-', 'LineWidth', 1.5);
                lg{end+1} = 'QZSS pooled 5 days'; %#ok<AGROW>
            end
            xlabel('Receiver bias offset (TECU)');
            if pp == 1, ylabel('Normalized cost (sum of VTEC std)'); end
            title(subt{pp});
            if pp == 2 && ~isempty(lg)
                legend(lg, 'Location', 'best');
            end
        end
        saveas(gcf, fullfile(outDir, 'dcb_costcurves.png'));
        close;
        fprintf('cost-curve diagnostic saved.\n');
    catch ME
        fprintf('cost-curve diagnostic skipped: %s\n', ME.message);
    end

    % ---------------- report + figure ----------------
    fid = fopen(fullfile(outDir, 'dcb_welch.txt'), 'w');
    fprintf(fid, 'RECEIVER-DCB / VTEC HYPOTHESIS TESTS\n');
    fprintf(fid, 't via incomplete beta; permutation p (10k shuffles, rng(42))\n\n');
    fprintf(fid, '%-42s %4s %10s %18s %8s %8s %10s %10s %7s\n', ...
        'test', 'n', 'estimate', '95% CI', 't', 'df', 'p(t)', 'p(perm)', 'd');
    for k = 1:numel(R)
        fprintf(fid, '%-42s %4d %+10.3f [%7.3f, %7.3f] %8.3f %8.1f %10.4g %10.4g %7.3f\n', ...
            R(k).test, R(k).n, R(k).est, R(k).ci(1), R(k).ci(2), ...
            R(k).t, R(k).df, R(k).p, R(k).pPerm, R(k).d);
        fprintf('%-42s n=%-3d est=%+8.3f p(t)=%.4g p(perm)=%.4g d=%.3f\n', ...
            R(k).test, R(k).n, R(k).est, R(k).p, R(k).pPerm, R(k).d);
    end
    fclose(fid);

    fig = figure('Name', 'welch', 'Position', [100 50 900, 220+110*max(numel(R),1)], 'Visible', 'off');
    hold on; grid on;
    for k = 1:numel(R)
        y = numel(R) - k + 1;
        plot(R(k).est, y, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        plot(R(k).ci, [y y], 'k-', 'LineWidth', 2);
    end
    xline(0, 'r--');
    set(gca, 'YTick', 1:numel(R), ...
        'YTickLabel', flip({R.test}), 'FontSize', 9);
    xlabel('Estimate [95% CI] (TECU)');
    title('Receiver-DCB / VTEC hypothesis tests');
    out = fullfile(outDir, 'dcb_welch.png');
    saveas(gcf, out);
    close;
    fprintf('welch saved: %s\n', out);
end

function r = pack(test, n, est, ci, t, df, p, pp, d)
    r = struct('test', test, 'n', n, 'est', est, 'ci', ci, ...
               't', t, 'df', df, 'p', p, 'pPerm', pp, 'd', d);
end

function [t, df, p, ci, dz, pp] = pairedTest(d)
% Paired t on difference vector d + permutation p + Cohen's dz.
    n = numel(d);
    m = mean(d); s = std(d);
    se = s / sqrt(n);
    t = m / se;
    df = n - 1;
    p = t2p(t, df);
    tc = tcrit(df);
    ci = [m - tc*se, m + tc*se];
    dz = m / s;
    pp = permP(d, 0);
end

function [t, df, p, ci, d, pp] = welchTest(x, y)
% Welch two-sample t + Satterthwaite df + Cohen's d + permutation p.
    nx = numel(x); ny = numel(y);
    mx = mean(x); my = mean(y);
    vx = var(x); vy = var(y);
    se = sqrt(vx/nx + vy/ny);
    t = (mx - my) / se;
    df = (vx/nx + vy/ny)^2 / ((vx/nx)^2/(nx-1) + (vy/ny)^2/(ny-1));
    p = t2p(t, df);
    tc = tcrit(df);
    ci = [(mx-my) - tc*se, (mx-my) + tc*se];
    sp = sqrt(((nx-1)*vx + (ny-1)*vy) / (nx+ny-2));
    d = (mx - my) / sp;
    pp = permP([x, y], nx);
end

function p = t2p(t, df)
% Two-sided t p-value via regularized incomplete beta (no toolbox).
    if ~isfinite(t) || df <= 0
        p = NaN; return;
    end
    x = df / (df + t^2);
    p = betainc(x, df/2, 0.5);
end

function tc = tcrit(df)
% Two-sided 95% t critical value. Uses built-in tinv when available,
% else Wilson-Hilferty/Cornish-Fisher expansion (accurate for df >= 3).
    try
        tc = tinv(0.975, df);
        return;
    catch
    end
    z = 1.959964;
    v = max(df, 3);
    tc = z * (1 + (z^2+1)/(4*v) + (5*z^4+16*z^2+3)/(96*v^2));
end

function pp = permP(varargin)
% Permutation p-value. pairedTest passes (d, 0): signs flipped.
% welchTest passes ([x y], nx): labels shuffled. 10k reps, rng(42).
    if nargin == 2 && isscalar(varargin{2}) && varargin{2} == 0
        d = varargin{1}(:);
        obs = abs(mean(d));
        s = rng; rng(42);
        cnt = 0; B = 10000;
        for b = 1:B
            sgn = sign(randn(numel(d), 1)); sgn(sgn == 0) = 1;
            if abs(mean(sgn .* d)) >= obs, cnt = cnt + 1; end
        end
        rng(s);
        pp = (cnt + 1) / (B + 1);
    else
        xy = varargin{1}(:); nx = varargin{2};
        obs = abs(mean(xy(1:nx)) - mean(xy(nx+1:end)));
        s = rng; rng(42);
        cnt = 0; B = 10000;
        for b = 1:B
            pr = xy(randperm(numel(xy)));
            if abs(mean(pr(1:nx)) - mean(pr(nx+1:end))) >= obs, cnt = cnt + 1; end
        end
        rng(s);
        pp = (cnt + 1) / (B + 1);
    end
end
