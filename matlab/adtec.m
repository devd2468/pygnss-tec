function d = adtec(x, winSec)
% ADTEC  Detrended TEC: subtract centered moving mean (NaN-aware).
% x: vector (per-second or per-epoch series), winSec: window in samples.
    if nargin < 2 || isempty(winSec)
        winSec = 3600;
    end
    d = nan(size(x));
    ok = ~isnan(x);
    if sum(ok) < 3
        return;
    end
    try
        tr = movmean(x, winSec, 'omitnan');
    catch
        % fallback for older MATLAB: manual moving mean ignoring NaN
        tr = nan(size(x));
        xi = find(ok);
        for k = 1:numel(xi)
            i = xi(k);
            j = xi(max(1, k-winSec):min(numel(xi), k+winSec));
            tr(i) = mean(x(j));
        end
    end
    d(ok) = x(ok) - tr(ok);
end
