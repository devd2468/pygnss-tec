function P = getIPP(dayRoot, station4, sysChar)
% GETIPP  Cached thin-shell IPPs in the day's Results folder.
% Computes via computeIPP on first call, loads thereafter.
% P fields: lat/lon/elev/azim/vtec/roti (86400 x NPRN), station/sys/...
    if nargin < 3 || isempty(sysChar), sysChar = 'G'; end
    resDir = fullfile(dayRoot, 'Results');
    f = fullfile(resDir, sprintf('IPP_%s_%s.mat', station4, sysChar));
    if exist(f, 'file')
        try
            P = load(f);
            if isfield(P, 'lat'), return; end
        catch
        end
    end
    P = computeIPP(dayRoot, station4, sysChar, resDir);
end
