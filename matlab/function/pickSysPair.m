function [PT, F, best] = pickSysPair(obsType, obsData, sys)
% PICKSYSPAIR  Best-validity dual-frequency pair for Galileo/BeiDou/QZSS.
% Tries all code-variant combinations (X/Q/L SVID flavors differ per
% receiver) and keeps the combo with the most jointly-valid observations.
% Output PT: struct(c1,c2,l1,l2); F: struct(f1,f2,k,lam1,lam2); best: count.
    PT = struct('c1', '', 'c2', '', 'l1', '', 'l2', '');
    F = struct('f1', NaN, 'f2', NaN, 'k', NaN, 'lam1', NaN, 'lam2', NaN);
    best = 0;
    c = 299792458; A = 40.3;
    switch sys
        case 'E'  % Galileo E1 + E5a (fallback E1 + E5b)
            groups = {
                1575.42e6, {'C1C','C1X'}, {'L1C','L1X'}, ...
                    1176.45e6, {'C5X','C5Q'}, {'L5X','L5Q'};
                1575.42e6, {'C1C','C1X'}, {'L1C','L1X'}, ...
                    1207.14e6, {'C7X','C7Q'}, {'L7X','L7Q'};
                };
        case 'C'  % BeiDou B1I + B2I (fallback B1C + B2a)
            groups = {
                1561.098e6, {'C2I'}, {'L2I'}, ...
                    1207.14e6, {'C7I'}, {'L7I'};
                1575.42e6, {'C1X'}, {'L1X'}, ...
                    1176.45e6, {'C5X'}, {'L5X'};
                };
        case 'J'  % QZSS L1 + L2 (fallback L1 + L5)
            groups = {
                1575.42e6, {'C1C','C1X'}, {'L1C','L1X'}, ...
                    1227.60e6, {'C2X','C2L'}, {'L2X','L2L'};
                1575.42e6, {'C1C','C1X'}, {'L1C','L1X'}, ...
                    1176.45e6, {'C5X','C5Q'}, {'L5X','L5Q'};
                };
        otherwise
            return;
    end
    ng = size(groups, 1);
    for g = 1:ng
        f1 = groups{g, 1}; C1 = groups{g, 2}; L1 = groups{g, 3};
        f2 = groups{g, 4}; C2 = groups{g, 5}; L2 = groups{g, 6};
        for a = 1:numel(C1)
            for b = 1:numel(L1)
                for d = 1:numel(C2)
                    for e = 1:numel(L2)
                        need = {C1{a}, C2{d}, L1{b}, L2{e}};
                        if ~all(ismember(need, obsType)), continue; end
                        A1 = obsData(:, ismember(obsType, need{1}));
                        B2 = obsData(:, ismember(obsType, need{2}));
                        P1 = obsData(:, ismember(obsType, need{3}));
                        P2 = obsData(:, ismember(obsType, need{4}));
                        nv = sum(~isnan(A1) & ~isnan(B2) & ~isnan(P1) & ~isnan(P2));
                        if nv > best
                            best = nv;
                            PT = struct('c1', need{1}, 'c2', need{2}, ...
                                        'l1', need{3}, 'l2', need{4});
                            k = f1^2*f2^2/(A*(f1^2-f2^2)*1e16);
                            F = struct('f1', f1, 'f2', f2, 'k', k, ...
                                       'lam1', c/f1, 'lam2', c/f2);
                        end
                    end
                end
            end
        end
    end
end
