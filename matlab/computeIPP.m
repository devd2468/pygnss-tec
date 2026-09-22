function P = computeIPP(dayRoot, station4, sysChar, outDir)
% COMPUTEIPP  True ionospheric pierce points (350 km thin shell) per
% epoch per satellite, recomputed from saved obs/nav products.
% Needs no reprocessing: reads MultiGNSS_<ST>[_<SYS>]_*.mat (VTEC/ROTI)
% plus the day folder's RINEX MO/MN files (same readers as the pipeline).
%
% Output IPP_<tag>.mat with fields lat/lon/elev/azim/vtec/roti
% (86400 x NPRN sparse matrices) + meta. NPRN=32 for G, 64 otherwise.
% Elevation is recomputed independently and cross-checked against the
% saved prm.elevation (logged; expect <0.1 deg).
    if nargin < 4 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if nargin < 3 || isempty(sysChar), sysChar = 'G'; end

    P = [];
    % ---- load TEC product ----
    if sysChar == 'G'
        pat = sprintf('MultiGNSS_%s_2*.mat', station4);
    else
        pat = sprintf('MultiGNSS_%s_%s_*.mat', station4, sysChar);
    end
    d = dir(fullfile(dayRoot, 'Results', pat));
    if isempty(d)
        % fallback: any matching file
        d = dir(fullfile(dayRoot, 'Results', sprintf('MultiGNSS_%s*.mat', station4)));
        keep = false(numel(d), 1);
        for k = 1:numel(d)
            hasSys = ~isempty(regexp(d(k).name, sprintf('_%s_\\d{4}_', sysChar), 'once'));
            isGPS = isempty(regexp(d(k).name, '_[ERCJ]_\d{4}_', 'once'));
            keep(k) = (sysChar == 'G' && isGPS) || (sysChar ~= 'G' && hasSys);
        end
        d = d(keep);
    end
    if isempty(d)
        fprintf('computeIPP: no product for %s/%s in %s\n', station4, sysChar, dayRoot);
        return;
    end
    S = load(fullfile(d(1).folder, d(1).name));
    V = S.TEC.vertical;
    R = [];
    if isfield(S, 'ROTI'), R = S.ROTI; end

    % ---- locate RINEX files ----
    mo = findMO(dayRoot, station4);
    mn = findNavFor(dayRoot, mo);
    if isempty(mo) || isempty(mn)
        fprintf('computeIPP: RINEX pair not found for %s in %s\n', station4, dayRoot);
        return;
    end

    % ---- read obs/nav per system ----
    if sysChar == 'G'
        [obs, nav] = readrinex304(mo, mn, [dayRoot filesep]);
        oType = obs.type; oData = obs.data;
        oIdx = obs.index; oEp = obs.epoch;
        eph = nav.eph; ephIdx = nav.index;
        NPRN = 32;
        rec = obs.rcvpos(:)';
    else
        O = rinex3obs_sys(fullfile(dayRoot, mo), sysChar);
        N = rinex3nav_ecj(fullfile(dayRoot, mn), sysChar);
        oType = O.Type.(sysChar); oData = O.Data.(sysChar)';
        oIdx = O.PRN.(sysChar)'; oEp = O.Ep.(sysChar)';
        eph = N.Eph.(sysChar)'; ephIdx = N.PRN.(sysChar)';
        NPRN = 64;
        rec = O.rcvpos(:)';
    end
    if isempty(rec) || any(~isfinite(rec))
        if isfield(S, 'obs') && isfield(S.obs, 'rcvpos')
            rec = S.obs.rcvpos(:)';
        else
            error('computeIPP:norec', 'No receiver position available.');
        end
    end
    ll0 = ecef2lla(rec);
    lat0 = ll0(1); lon0 = ll0(2);

    % ENU basis at receiver (same convention as PositionA2B)
    Ev = [-sind(lon0), cosd(lon0), 0];
    Nv = [-sind(lat0)*cosd(lon0), -sind(lat0)*sind(lon0), cosd(lat0)];
    Uv = [cosd(lat0)*cosd(lon0), cosd(lat0)*sind(lon0), sind(lat0)];

    % code column for travel-time correction (any C* type)
    ci = find(cellfun(@(t) ~isempty(t) && t(1) == 'C', oType), 1);
    if isempty(ci), ci = 1; end

    nR = size(V, 1);
    lat = nan(nR, NPRN); lon = nan(nR, NPRN);
    elv = nan(nR, NPRN); azi = nan(nR, NPRN);
    prns = unique(oIdx);
    prns = prns(prns >= 1 & prns <= NPRN);
    maxElErr = 0;
    for i = 1:numel(prns)
        p = prns(i);
        rows = find(oIdx == p);
        tm0 = round(oEp(rows));
        keep = tm0 >= 0 & tm0 < nR;
        rows = rows(keep);
        if isempty(rows), continue; end
        vv = V(tm0(keep)+1, p);
        okV = isfinite(vv);
        if sum(okV) < 3, continue; end
        rows = rows(okV); tm = tm0(keep); tm = tm(okV);
        ps = oData(rows, ci);
        pv = ps(isfinite(ps));
        if isempty(pv), ps(:) = 22e6; else, ps(~isfinite(ps)) = median(pv); end
        try
            [sp, ~] = satpos_xyz_sbias(tm, p, eph, ephIdx, ps);
        catch
            continue;
        end
        rv = sp - rec;
        nr = sqrt(sum(rv.^2, 2));
        rh = rv ./ nr;
        e = asind(rh * Uv');
        a = atan2d(rh * Ev', rh * Nv');
        % thin-shell IPP, h = 350 km
        [la, lo] = ippFromElAz(lat0, lon0, e, a, 350);
        ii = tm + 1;
        lat(ii, p) = la; lon(ii, p) = lo;
        elv(ii, p) = e; azi(ii, p) = a;
    end
    % elevation cross-check vs saved product
    if isfield(S, 'prm') && isfield(S.prm, 'elevation')
        E0 = S.prm.elevation;
        m = min(size(E0, 1), nR);
        w = min(size(E0, 2), NPRN);
        dd = abs(elv(1:m, 1:w) - E0(1:m, 1:w));
        dd = dd(isfinite(dd));
        if ~isempty(dd), maxElErr = max(dd); end
    end
    fprintf('computeIPP %s/%s: IPPs=%d, elev-xcheck max err=%.3f deg\n', ...
        station4, sysChar, sum(isfinite(lat(:))), maxElErr);

    P = struct('lat', lat, 'lon', lon, 'elev', elv, 'azim', azi, ...
        'vtec', V(:, 1:min(NPRN, size(V, 2))), 'roti', [], ...
        'station', station4, 'sys', sysChar, 'reflat', lat0, 'reflon', lon0);
    if ~isempty(R)
        P.roti = R(:, 1:min(NPRN, size(R, 2)));
    end
    out = fullfile(outDir, sprintf('IPP_%s_%s.mat', station4, sysChar));
    save(out, '-struct', 'P');
    fprintf('computeIPP saved: %s\n', out);
end

function [la, lo] = ippFromElAz(lat0, lon0, E, A, h)
% Thin-shell IPP via great-circle destination formula.
% lat0/lon0 deg, E elevation deg, A azimuth deg clockwise from north.
    R = 6371.0;
    lat0 = deg2rad(lat0); lon0 = deg2rad(lon0);
    E = deg2rad(E); A = deg2rad(A);
    psi = pi/2 - E - asin(R * cos(E) / (R + h));  % central angle
    la = asin(sin(lat0) .* cos(psi) + cos(lat0) .* sin(psi) .* cos(A));
    lo = lon0 + atan2(sin(A) .* sin(psi) .* cos(lat0), ...
                      cos(psi) - sin(lat0) .* sin(la));
    la = rad2deg(la); lo = rad2deg(lo);
    lo = mod(lo + 540, 360) - 180;
end

function mo = findMO(dataRoot, station4)
    mo = '';
    d = dir(fullfile(dataRoot, '*.rnx'));
    for k = 1:numel(d)
        nm = d(k).name;
        if strncmp(nm, station4, 4) && isempty(regexp(nm, '_[MEGRSI]N\.rnx$', 'once')) ...
                && ~isempty(regexp(nm, '(_MO\.rnx|\.\d\do)$', 'once'))
            mo = nm;
            return;
        end
    end
end

function mn = findNavFor(dataRoot, moName)
    mn = '';
    if isempty(moName), return; end
    t = regexprep(moName, '_\d+S_MO\.rnx$', '_MN.rnx', 'ignorecase');
    if ~strcmp(t, moName) && exist(fullfile(dataRoot, t), 'file')
        mn = t; return;
    end
    t = regexprep(moName, '_MO\.rnx$', '_MN.rnx', 'ignorecase');
    if ~strcmp(t, moName) && exist(fullfile(dataRoot, t), 'file')
        mn = t; return;
    end
    t = regexprep(moName, '\.([0-9A-Za-z]{2})o$', '.$1n');
    if ~strcmp(t, moName) && exist(fullfile(dataRoot, t), 'file')
        mn = t; return;
    end
end
