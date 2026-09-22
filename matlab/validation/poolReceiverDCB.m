function [rcv, costCurve, brGrid, pegged] = poolReceiverDCB(dayRoots, station4, sysChar)
% POOLRECEIVERDCB  Multi-day joint Ma-Maruyama receiver-bias estimation.
% Pools withrcvbias + slant factor across days for one station/system,
% then runs a single minimum-spread search. Cures single-day sparsity
% (e.g. QZSS 3-4 sats pegging the search bounds) without any circularity:
% the pooled estimate uses only that system's own data.
%
% The search runs directly in TECU-offset space (equivalent to the ns
% search in rcv_bias_sys up to the monotone ns->TECU map), so no per-pair
% frequencies are needed and the result compares directly with DCB.rcv.
%
% Outputs:
%   rcv       receiver bias (TECU), NaN if insufficient data
%   costCurve sum-of-std vs brGrid (for diagnostic plots)
%   brGrid    searched bias grid (TECU)
%   pegged    true if minimum sits at a search bound (unreliable estimate)
    rcv = NaN; costCurve = []; brGrid = []; pegged = false;

    ST = []; SF = [];
    for d = 1:numel(dayRoots)
        if strcmp(sysChar, 'G')
            pat = sprintf('MultiGNSS_%s_2*.mat', station4);
        else
            pat = sprintf('MultiGNSS_%s_%s_*.mat', station4, sysChar);
        end
        dd = dir(fullfile(dayRoots{d}, 'Results', pat));
        % exclude cross-system files for GPS pattern
        keep = true(numel(dd), 1);
        if strcmp(sysChar, 'G')
            for k = 1:numel(dd)
                if ~isempty(regexp(dd(k).name, '^MultiGNSS_[A-Z0-9]{4}_[ERCJ]_', 'once'))
                    keep(k) = false;
                end
            end
        end
        dd = dd(keep);
        if isempty(dd), continue; end
        try
            S = load(fullfile(dd(1).folder, dd(1).name));
            W = S.TEC.withrcvbias;
            E = S.prm.elevation;
            Re = 6371.009e3; h = 350e3;
            sf = sqrt(1 - (Re * cosd(E) ./ (Re + h)).^2);
            ST = cat(1, ST, W); %#ok<AGROW>
            SF = cat(1, SF, sf); %#ok<AGROW>
        catch
        end
    end
    % downsample like rcv_bias_sys (1:30) to bound cost
    ST = ST(1:30:end, :);
    SF = SF(1:30:end, :);
    if size(ST, 1) < 100 || sum(isfinite(ST(:))) < 1000
        return;
    end
    brGrid = -90:1.5:90;  % TECU offsets; coarse grid for diagnostics
    costCurve = nan(size(brGrid));
    for k = 1:numel(brGrid)
        V = (ST - brGrid(k)) .* SF;
        costCurve(k) = nansum(nanstd(V'));
    end
    % refine around coarse minimum (two refinement passes like rcv_bias_sys)
    [~, ix] = nanmin(costCurve);
    if isempty(ix) || all(isnan(costCurve))
        return;
    end
    lo = max(brGrid(1), brGrid(ix) - 3);
    hi = min(brGrid(end), brGrid(ix) + 3);
    g = brGrid; cc = costCurve;
    for pass = 1:2
        g = linspace(lo, hi, 41);
        cc = nan(size(g));
        for k = 1:numel(g)
            V = (ST - g(k)) .* SF;
            cc(k) = nansum(nanstd(V'));
        end
        [~, ix2] = nanmin(cc);
        lo = max(brGrid(1), g(max(ix2-1, 1)));
        hi = min(brGrid(end), g(min(ix2+1, numel(g))));
    end
    [~, ix3] = nanmin(cc);
    rcv = g(ix3);
    pegged = (ix3 == 1) || (ix3 == numel(g));
end
