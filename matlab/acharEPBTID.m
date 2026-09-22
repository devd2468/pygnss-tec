function evts = acharEPBTID(dayRoots, stations, outDir)
% ACHAREPBTID  EPB/TID characterization for the results section.
% Detects ionospheric irregularity events per station-day-PRN, estimates
% peak ROTI, DTEC amplitude/depletion, dominant wave period (FFT), and
% classifies EPB vs TID candidates using local-time + morphology rules.
%
% Outputs:
%   epb_tid_table.csv/.txt  event catalog (results-section ready)
%   acharEvent_<ST>_<date>_PRN<nn>.png  per-event zoom figures
%   acharTimeline.png        all-events timeline across days
%
% Rules (transparent, adjustable below):
%   EPB candidate : local time 17.5-24h IST AND (peakROTI>=0.5 | depth<=-5)
%   TID candidate : spectral peak 15-90 min with prominence>=2
%   else          : irregularity / weak fluctuation by peakROTI level
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    TH_ROTI = 0.3;    % detection floor (TECU/min)
    TH_PEAK = 0.5;    % strong-event level
    MIN_DUR = 10;     % minutes
    MERGE_GAP = 5;    % minutes
    IST = 5.5;        % UTC+5:30

    evts = struct('day', {}, 'station', {}, 'prn', {}, ...
        't0utc', {}, 't1utc', {}, 't0lt', {}, 't1lt', {}, ...
        'peakROTI', {}, 'tPeakUTC', {}, 'dtecPTP', {}, 'depth', {}, ...
        'domPeriod', {}, 'prom', {}, 'nSup', {}, 'class', {});

    for d = 1:numel(dayRoots)
        dayTag = sprintf('day%02d', d);
        try
            dn = datetime(regexp(dayRoots{d}, '\d{4}-\d{2}-\d{2}', 'match', 'once'));
        catch
            dn = NaT;
        end
        for s = 1:numel(stations)
            st = stations{s};
            S = aloadTEC(dayRoots{d}, st);
            if isempty(S) || ~isfield(S, 'TEC') || ~isfield(S, 'ROTI')
                continue;
            end
            V = S.TEC.vertical;
            R = S.ROTI;
            nPRN = size(V, 2);
            for p = 1:nPRN
                r = R(:, p);
                if sum(isfinite(r)) < 60, continue; end
                % contiguous ROTI>=TH segments (per-second rows)
                m = isfinite(r) & r >= TH_ROTI;
                dm = diff([false; m; false]);
                stt = find(dm == 1); enn = find(dm == -1) - 1;
                % merge small gaps
                k = 1;
                while k < numel(stt)
                    if (stt(k+1) - enn(k)) / 60 <= MERGE_GAP
                        enn(k) = enn(k+1);
                        stt(k+1) = []; enn(k+1) = [];
                    else
                        k = k + 1;
                    end
                end
                for k = 1:numel(stt)
                    dur = (enn(k) - stt(k) + 1) / 60;
                    if dur < MIN_DUR, continue; end
                    seg = stt(k):enn(k);
                    [pk, ip] = nanmax(r(seg));
                    if ~isfinite(pk), continue; end
                    vs = V(seg, p);
                    dd = adtec(V(:, p), 3600);
                    dw = dd(seg);
                    ptp = nanmax(dw) - nanmin(dw);
                    if ~isfinite(ptp), ptp = NaN; end
                    dep = nanmin(dw);
                    % dominant period via FFT on 30 s resample
                    [per, prom] = domPeriod(tHvec(seg), dw);
                    t0u = (stt(k)-1)/3600; t1u = (enn(k)-1)/3600;
                    tpk = (seg(ip)-1)/3600;
                    % supporting sats (other PRNs active in same window)
                    nsup = 0;
                    for q = 1:nPRN
                        if q == p, continue; end
                        rq = R(seg, q);
                        if sum(isfinite(rq) & rq >= TH_ROTI) >= 5, nsup = nsup + 1; end
                    end
                    cls = classify(pk, dep, tpk + IST, per, prom);
                    evts(end+1) = struct('day', dayTag, 'station', st, ...
                        'prn', p, 't0utc', t0u, 't1utc', t1u, ...
                        't0lt', t0u + IST, 't1lt', t1u + IST, ...
                        'peakROTI', pk, 'tPeakUTC', tpk, ...
                        'dtecPTP', ptp, 'depth', dep, ...
                        'domPeriod', per, 'prom', prom, ...
                        'nSup', nsup, 'class', cls); %#ok<AGROW>
                end
            end
        end
    end

    if isempty(evts)
        fprintf('acharEPBTID: no events detected.\n');
        return;
    end

    % keep top events per station-day by peak ROTI (cap figure count)
    try
        keys = strcat({evts.day}, '_', {evts.station});
        [~, ~, ic] = unique(keys);
        keep = false(size(evts));
        for k = 1:max(ic)
            ix = find(ic == k);
            [~, o] = sort([evts(ix).peakROTI], 'descend');
            keep(ix(o(1:min(3, numel(o))))) = true;
        end
    catch
        keep = true(size(evts));
    end

    % ---- write catalog ----
    csv = fullfile(outDir, 'epb_tid_table.csv');
    fid = fopen(csv, 'w');
    fprintf(fid, 'day,station,prn,t0utc,t1utc,t0lt,t1lt,peakROTI,tPeakUTC,dtecPTP,depth,domPeriod_min,prominence,nSupport,class\n');
    for k = 1:numel(evts)
        e = evts(k);
        fprintf(fid, '%s,%s,G%02d,%.2f,%.2f,%.2f,%.2f,%.3f,%.2f,%.2f,%.2f,%.1f,%.2f,%d,%s\n', ...
            e.day, e.station, e.prn, e.t0utc, e.t1utc, e.t0lt, e.t1lt, ...
            e.peakROTI, e.tPeakUTC, e.dtecPTP, e.depth, e.domPeriod, ...
            e.prom, e.nSup, e.class);
    end
    fclose(fid);
    % txt summary
    txt = fullfile(outDir, 'epb_tid_table.txt');
    fid = fopen(txt, 'w');
    fprintf(fid, 'EPB/TID EVENT CATALOG\n');
    fprintf(fid, 'Rules: EPB=LT 17.5-24 IST & (peakROTI>=%.1f | depth<=-5); TID=spectral 15-90 min & prom>=2\n', TH_PEAK);
    fprintf(fid, 'Total candidate windows: %d | plotted (top-3/station-day): %d\n\n', numel(evts), sum(keep));
    fprintf(fid, '%-6s %-7s %-4s %-11s %-9s %-8s %-8s %-6s %s\n', ...
        'day', 'station', 'prn', 'UTC', 'LT(IST)', 'pkROTI', 'ptp', 'per', 'class');
    for k = 1:numel(evts)
        e = evts(k);
        fprintf(fid, '%-6s %-7s G%02d %05.2f-%05.2f %05.2f-%05.2f %-8.3f %-8.2f %-6.1f %s%s\n', ...
            e.day, e.station, e.prn, e.t0utc, e.t1utc, e.t0lt, e.t1lt, ...
            e.peakROTI, e.dtecPTP, e.domPeriod, e.class, ...
            ternary(keep(k), ' *', ''));
    end
    fprintf(fid, '\n* = plotted figure\n');
    fclose(fid);
    fprintf('acharEPBTID: %d events -> %s\n', numel(evts), csv);

    % ---- per-event figures (kept subset) ----
    ki = find(keep);
    for k = 1:numel(ki)
        try
            plotEventFig(dayRoots, evts(ki(k)), outDir);
        catch ME
            fprintf('  event fig failed: %s\n', ME.message);
        end
    end

    % ---- timeline ----
    try
        plotTimeline(evts, outDir);
    catch ME
        fprintf('  timeline failed: %s\n', ME.message);
    end
end

function cls = classify(pk, dep, ltPk, per, prom)
    if ltPk >= 17.5 && ltPk <= 24.5 && (pk >= 0.5 || (isfinite(dep) && dep <= -5))
        cls = 'EPB candidate';
    elseif isfinite(per) && per >= 15 && per <= 90 && prom >= 2
        cls = 'TID candidate';
    elseif pk >= 0.5
        cls = 'Irregularity';
    elseif pk >= 0.3
        cls = 'Weak fluctuation';
    else
        cls = 'Quiet';
    end
end

function [per, prom] = domPeriod(tH, dw)
% Dominant wave period (min) in 10-120 min band via FFT on 30 s grid.
    per = NaN; prom = NaN;
    try
        ok = isfinite(dw) & isfinite(tH(:));
        if sum(ok) < 20, return; end
        t = tH(ok); y = dw(ok);
        y = y - median(y);
        % 30 s uniform grid over span
        tg = (t(1):30/3600:t(end))';
        if numel(tg) < 20, return; end
        yg = interp1(t, y, tg, 'linear', NaN);
        okg = isfinite(yg);
        if sum(okg) < 20, return; end
        yg(~okg) = 0;  % zero-fill short gaps after check
        Y = abs(fft(yg - mean(yg)));
        n = numel(Y);
        P = (n*30/60) ./ (1:floor(n/2));  % periods in min
        Pw = Y(1:floor(n/2));
        sel = P >= 10 & P <= 120;
        if sum(sel) < 3, return; end
        [mx, ix] = max(Pw(sel));
        ii = find(sel); ii = ii(ix);
        per = P(ii);
        prom = mx / (median(Pw(sel)) + eps);
    catch
    end
end

function t = tHvec(seg)
    t = (seg - 1) / 3600;
end

function o = ternary(c, a, b)
    if c, o = a; else, o = b; end
end

function plotEventFig(dayRoots, e, outDir)
    dnum = sscanf(e.day, 'day%d');
    S = aloadTEC(dayRoots{dnum}, e.station);
    V = S.TEC.vertical(:, e.prn);
    R = S.ROTI(:, e.prn);
    tH = ((0:numel(V)-1) / 3600)';  % column: matches data orientation
    pad = 0.75;
    w0 = max(tH(1), e.t0utc - pad); w1 = min(tH(end), e.t1utc + pad);
    inW = tH >= w0 & tH <= w1;
    d = adtec(V, 3600);
    d = d(:); V = V(:); R = R(:);

    fig = figure('Name', 'event', 'Position', [100 50 1000 900], 'Visible', 'off');
    try
        ax1 = subplot(3, 1, 1); hold on; grid on;
        ok = isfinite(V) & inW;
        plot(tH(ok), V(ok), 'k-', 'LineWidth', 1.2);
        ylabel('VTEC (TECU)');
        title(sprintf('%s %s G%02d  %s  pkROTI=%.2f  ptp=%.2f TECU  T~%.0f min', ...
            e.day, e.station, e.prn, e.class, e.peakROTI, e.dtecPTP, e.domPeriod));
        ax2 = subplot(3, 1, 2); hold on; grid on;
        ok = isfinite(d) & inW;
        plot(tH(ok), d(ok), 'b-', 'LineWidth', 1.1);
        yline(0, 'k--');
        ylabel('DTEC (TECU)');
        ax3 = subplot(3, 1, 3); hold on; grid on;
        ok = isfinite(R) & inW;
        plot(tH(ok), R(ok), 'r.-', 'MarkerSize', 5);
        yline(0.2, 'b--'); yline(0.5, 'r--');
        ylabel('ROTI (TECU/min)');
        xlabel(sprintf('Universal Time (hour)  [IST = UT+5:30; event %.2f-%.2f UT]', e.t0utc, e.t1utc));
        for a = [ax1, ax2, ax3]
            xlim(a, [w0, w1]);
        end
        out = fullfile(outDir, sprintf('acharEvent_%s_%s_G%02d.png', e.day, e.station, e.prn));
        saveas(gcf, out);
        fprintf('  event fig: %s\n', out);
    catch ME
        fprintf('  event fig failed: %s\n', ME.message);
    end
    try, close(fig); catch, end
end

function plotTimeline(evts, outDir)
    days = unique({evts.day});
    sts = unique({evts.station});
    fig = figure('Name', 'timeline', 'Position', [100 50 1200 200+90*numel(sts)], 'Visible', 'off');
    hold on; grid on;
    cols = containers.Map({'EPB candidate', 'TID candidate', 'Irregularity', 'Weak fluctuation', 'Quiet'}, ...
        {[0.85 0.1 0.1], [0.1 0.4 0.85], [1 0.55 0], [0.4 0.7 0.4], [0.7 0.7 0.7]});
    for s = 1:numel(sts)
        for k = 1:numel(evts)
            e = evts(k);
            if ~strcmp(e.station, sts{s}), continue; end
            dn = sscanf(e.day, 'day%d');
            x0 = (dn-1)*24 + e.t0utc; x1 = (dn-1)*24 + e.t1utc;
            y0 = s - 0.35; y1 = s + 0.35;
            try, cc = cols(e.class); catch, cc = [0.5 0.5 0.5]; end
            patch([x0 x1 x1 x0], [y0 y0 y1 y1], cc, 'EdgeColor', 'k');
        end
    end
    set(gca, 'YTick', 1:numel(sts), 'YTickLabel', sts);
    xlabel('Hours past 01-May 00:00 UT');
    ylabel('Station');
    title('Irregularity event timeline (color = class, width = duration)');
    xlim([0, numel(days)*24]);
    out = fullfile(outDir, 'acharTimeline.png');
    saveas(gcf, out);
    close;
    fprintf('  timeline: %s\n', out);
end
