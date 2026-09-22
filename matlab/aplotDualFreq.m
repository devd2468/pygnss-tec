function aplotDualFreq(dayRoot, station4, outDir)
% APLOTDUALFREQ  Dual-frequency ROTI scatter (BDS-style figure, GPS data).
% Compares GPS L1-L2 vs L1-L5 carrier ROTI, recomputed from saved obs.
% Colors: L1-L2 red, L1-L5 cyan (cf. L1P-L5P / L2I-L6I panels).
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    addpathSafe();

    S = aloadTEC(dayRoot, station4);
    if isempty(S) || ~isfield(S, 'obs')
        fprintf('aplotDualFreq: no obs for %s in %s\n', station4, dayRoot);
        return;
    end
    obs = S.obs;
    c = 299792458; A = 40.3;
    f1 = 1575.42e6; f2 = 1227.60e6; f5 = 1176.45e6;
    k12 = f1^2*f2^2/(A*(f1^2-f2^2)*1e16);
    k15 = f1^2*f5^2/(A*(f1^2-f5^2)*1e16);
    lam1 = c/f1; lam2 = c/f2; lam5 = c/f5;

    % pick best-validity L1/L2 pair (same rule as engine), plus L5 if present
    [tC1, tL1, tC2, tL2] = pickBestL12(obs);
    tC5 = pickFirst(obs.type, {'C5X', 'C5Q'});
    tL5 = pickFirst(obs.type, {'L5X', 'L5Q'});
    if isempty(tC1) || isempty(tL1) || isempty(tC2) || isempty(tL2)
        fprintf('aplotDualFreq: L1/L2 pair unavailable for %s\n', station4);
        return;
    end
    hasL5 = ~isempty(tC5) && ~isempty(tL5);

    sats = unique(obs.index);
    sats = sats(sats >= 1 & sats <= 32);
    ST12 = nan(86400, 32); ST15 = nan(86400, 32);
    for i = 1:numel(sats)
        p = sats(i);
        rows = find(obs.index == p);
        tm = round(obs.epoch(rows));
        ok = tm >= 0 & tm < 86400;
        rows = rows(ok); tm = tm(ok);
        if isempty(rows), continue; end
        C1 = obs.data(rows, ismember(obs.type, tC1));
        P2 = obs.data(rows, ismember(obs.type, tC2));
        L1 = obs.data(rows, ismember(obs.type, tL1));
        L2 = obs.data(rows, ismember(obs.type, tL2));
        good = ~isnan(L1) & ~isnan(L2);
        ST12(tm(good)+1, p) = k12 * (lam1*L1(good) - lam2*L2(good));
        if hasL5
            C5 = obs.data(rows, ismember(obs.type, tC5));
            L5 = obs.data(rows, ismember(obs.type, tL5));
            good5 = ~isnan(L1) & ~isnan(L5) & ~isnan(C5);
            ST15(tm(good5)+1, p) = k15 * (lam1*L1(good5) - lam5*L5(good5));
            %#ok<NASGU> % C-fields unused: phase-only ROTI
        end
    end
    R12 = roticalculation(ST12, sats');
    R15 = hasL5 * roticalculation(ST15, sats'); %#ok<NASGU>
    if ~hasL5, R15 = nan(size(R12)); end

    fig = figure('Name', sprintf('Dual-freq %s', station4), ...
        'Position', [100 50 1000 700], 'Visible', 'off');
    hold on; grid on;
    tH = (0:86399) / 3600;
    % decimate for scatter speed
    h12 = plotScatterDec(tH, R12, [1 0 0]);
    hleg = {'L1-L2'};
    if hasL5
        h15 = plotScatterDec(tH, R15, [0 0.8 0.9]);
        hleg{end+1} = 'L1-L5'; %#ok<AGROW>
    end
    xlabel('Universal Time (hour)'); ylabel('ROTI (TECU/min)');
    title(sprintf('%s  GPS dual-frequency ROTI  (%s)', station4, ...
        datestr(obsDate(obs), 'dd-mmm-yyyy')));
    legend(hleg, 'Location', 'northeast');
    ylim([0, max(1, nanmax(R12(:)) * 1.1)]);
    xlim([0, 24]);
    out = fullfile(outDir, sprintf('aplotDualFreq_%s.png', station4));
    saveas(gcf, out);
    close;
    fprintf('aplotDualFreq saved: %s\n', out);
end

function h = plotScatterDec(tH, R, col)
    h = [];
    [nr, nc] = size(R);
    keep = 1:10:nr;  % decimate rows
    xx = []; yy = [];
    for p = 1:nc
        v = R(keep, p);
        ok = ~isnan(v);
        if any(ok)
            xx = [xx; tH(keep(ok))']; %#ok<AGROW>
            yy = [yy; v(ok)]; %#ok<AGROW>
        end
    end
    if ~isempty(xx)
        h = scatter(xx, yy, 6, col, 'filled', 'MarkerFaceAlpha', 0.5);
    end
end

function t = pickFirst(types, cands)
    t = '';
    for k = 1:numel(cands)
        if any(strcmp(types, cands{k})), t = cands{k}; return; end
    end
end

function [tC1, tL1, tC2, tL2] = pickBestL12(obs)
    tC1 = ''; tL1 = ''; tC2 = ''; tL2 = '';
    cands = {
        'C1C', 'C2W', 'L1C', 'L2W';
        'C1C', 'C2X', 'L1C', 'L2X';
        'C1X', 'C2X', 'L1X', 'L2X';
        'C1X', 'C2W', 'L1X', 'L2W';
        };
    best = 0;
    for k = 1:size(cands, 1)
        need = cands(k, :);
        if all(ismember(need, obs.type))
            C = obs.data(:, ismember(obs.type, need{1}));
            P = obs.data(:, ismember(obs.type, need{2}));
            A1 = obs.data(:, ismember(obs.type, need{3}));
            A2 = obs.data(:, ismember(obs.type, need{4}));
            nv = sum(~isnan(C) & ~isnan(P) & ~isnan(A1) & ~isnan(A2));
            if nv > best
                best = nv;
                tC1 = need{1}; tC2 = need{2}; tL1 = need{3}; tL2 = need{4};
            end
        end
    end
end

function d = obsDate(obs)
    try
        d = datetime(obs.date(1), obs.date(2), obs.date(3));
    catch
        d = datetime('today');
    end
end

function addpathSafe()
    try
        p = fileparts(mfilename('fullpath'));
        addpath(p);
        if exist(fullfile(p, 'function'), 'dir')
            addpath(fullfile(p, 'function'));
        end
    catch
    end
end
