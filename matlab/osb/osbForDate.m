function osb = osbForDate(MGEX_ROOT, ymd)
% OSBFORDATE  Cached readCODEosb for one calendar date (persistent cache).
% ymd: [year month day]. Returns [] if no matching OSB file found.
    persistent cache;  % struct array: key, osb
    if isempty(cache), cache = struct('key', {}, 'osb', {}); end
    osb = [];
    if isempty(MGEX_ROOT) || ~exist(MGEX_ROOT, 'dir')
        return;
    end
    dn = datenum(ymd);
    doy = round(dn - datenum(ymd(1), 1, 1)) + 1;
    key = sprintf('%04d%03d', ymd(1), doy);
    for k = 1:numel(cache)
        if strcmp(cache(k).key, key)
            osb = cache(k).osb;
            return;
        end
    end
    d = dir(fullfile(MGEX_ROOT, sprintf('COD0MGXFIN_%s*OSB.BIA', key)));
    if isempty(d)
        return;
    end
    try
        o = readCODEosb(fullfile(d(1).folder, d(1).name));
        cache(end+1) = struct('key', key, 'osb', o); %#ok<AGROW>
        osb = o;
    catch ME
        fprintf('  OSB read note: %s\n', ME.message);
    end
end
