function [tC1, tP2, tL1, tL2] = pickFreqPair(obs)
% PICKFREQPAIR  Pick the GPS L1/L2 code+phase pair with the most
% jointly-valid observations. Shared by the TEC engine and the OSB-DCB
% builder so both use identical signal selection.
% Candidates cover legacy (C2W/L2W) and modern (C2X/L2X) receivers.
cands = {
    'C1C', 'C2W', 'L1C', 'L2W';
    'C1C', 'C2X', 'L1C', 'L2X';
    'C1X', 'C2X', 'L1X', 'L2X';
    'C1X', 'C2W', 'L1X', 'L2W';
    };
best = 0; bi = 1;
for k = 1:size(cands, 1)
    need = cands(k, :);
    if all(ismember(need, obs.type))
        C = obs.data(:, ismember(obs.type, need{1}));
        P = obs.data(:, ismember(obs.type, need{2}));
        L1v = obs.data(:, ismember(obs.type, need{3}));
        L2v = obs.data(:, ismember(obs.type, need{4}));
        nv = sum(~isnan(C) & ~isnan(P) & ~isnan(L1v) & ~isnan(L2v));
        if nv > best, best = nv; bi = k; end
    end
end
tC1 = cands{bi, 1}; tP2 = cands{bi, 2};
tL1 = cands{bi, 3}; tL2 = cands{bi, 4};
end
