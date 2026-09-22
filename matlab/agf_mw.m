function [GF, MW, prns, tag] = agf_mw(obs)
% AGF_MW  Geometry-free phase (m) and Melbourne-Wuebbena (cycles) per PRN.
% Uses the best-available L1/L2 pair in obs (same priority as the engine).
% GF, MW are 86400x32 matrices aligned to Time+1 rows.
    GF = nan(86400, 32); MW = nan(86400, 32);
    prns = []; tag = '';
    c = 299792458;
    f1 = 1575.42e6; f2 = 1227.60e6;
    lam1 = c/f1; lam2 = c/f2;
    pairs = {
        'C1C', 'C2W', 'L1C', 'L2W';
        'C1C', 'C2X', 'L1C', 'L2X';
        'C1X', 'C2X', 'L1X', 'L2X';
        'C1X', 'C2W', 'L1X', 'L2W';
        };
    % best joint-validity pair (same rule as the TEC engine)
    bi = 0; best = 0;
    for k = 1:size(pairs, 1)
        need = pairs(k, :);
        if all(ismember(need, obs.type))
            C = obs.data(:, ismember(obs.type, need{1}));
            P = obs.data(:, ismember(obs.type, need{2}));
            A1 = obs.data(:, ismember(obs.type, need{3}));
            A2 = obs.data(:, ismember(obs.type, need{4}));
            nv = sum(~isnan(C) & ~isnan(P) & ~isnan(A1) & ~isnan(A2));
            if nv > best, best = nv; bi = k; end
        end
    end
    if bi == 0, return; end
    tC1 = pairs{bi,1}; tP2 = pairs{bi,2};
    tL1 = pairs{bi,3}; tL2 = pairs{bi,4};
    tag = sprintf('%s/%s+%s/%s', tC1, tL1, tP2, tL2);
    sats = unique(obs.index);
    sats = sats(sats >= 1 & sats <= 32);
    aw = (f1-f2)/(f1+f2);
    for i = 1:numel(sats)
        p = sats(i);
        rows = find(obs.index == p);
        tm = round(obs.epoch(rows));
        ok = tm >= 0 & tm < 86400;
        rows = rows(ok); tm = tm(ok);
        if isempty(rows), continue; end
        P1 = obs.data(rows, ismember(obs.type, tC1));
        P2 = obs.data(rows, ismember(obs.type, tP2));
        L1 = obs.data(rows, ismember(obs.type, tL1));
        L2 = obs.data(rows, ismember(obs.type, tL2));
        good = ~isnan(P1) & ~isnan(P2) & ~isnan(L1) & ~isnan(L2);
        if ~any(good), continue; end
        ii = tm(good)+1;
        GF(ii, p) = lam1*L1(good) - lam2*L2(good);
        MW(ii, p) = (L1(good)-L2(good)) - aw*(P1(good)/lam1 + P2(good)/lam2);
        prns(end+1) = p; %#ok<AGROW>
    end
end
