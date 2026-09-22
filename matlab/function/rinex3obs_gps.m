function Obs = rinex3obs_gps(obsFullPath)
% RINEX3OBS_GPS  Pure-MATLAB RINEX 3.x observation reader (GPS only).
% Replaces the crashing ReadOBSrinex304.mexw64 binary.
%
% Input : obsFullPath - full path to RINEX 3 observation file
% Output: Obs struct with fields mimicking the MEX output used by
%         readrinex304.m:
%           Obs.Com      - {} (unused)
%           Obs.XYZr     - 1x3 approx receiver position (m)
%           Obs.Type.G   - 1xM cell of GPS observation codes
%           Obs.Date.St  - 1x6 first-obs date [y m d h mi s]
%           Obs.Ep.G     - 1xN seconds-of-day per GPS observation row
%           Obs.Data.G   - MxN observation matrix (cols = rows)
%           Obs.PRN.G    - 1xN GPS PRN numbers
%           Obs.HeaderSys- constellations declared in header
%
% Only GPS ('G') satellite lines are stored; other systems are counted
% from the header for coverage logging.

Obs = struct('Com', {{}}, 'XYZr', [0 0 0], 'HeaderSys', {{}});

fid = fopen(obsFullPath, 'r');
if fid < 0, error('Cannot open observation file: %s', obsFullPath); end

% ---------- header ----------
types = struct(); order = {};
xyzr = [0 0 0]; firstObs = [];
inHeader = true;
pendingTypes = {}; pendingNeed = 0;
while inHeader
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if length(ln) < 60, ln = [ln repmat(' ', 1, 60-length(ln))]; end
    label = strtrim(ln(61:min(80, length(ln))));
    % continuation of SYS types from previous line
    if pendingNeed > 0
        toks = strsplit(strtrim(ln(1:60)));
        toks = toks(~cellfun('isempty', toks));
        pendingTypes = [pendingTypes, toks]; %#ok<AGROW>
        if numel(pendingTypes) >= pendingNeed
            types.(pendingSys) = pendingTypes(1:pendingNeed);
            order{end+1} = pendingSys; %#ok<AGROW>
            pendingNeed = 0; pendingTypes = {};
        end
        continue;
    end
    if ~isempty(strfind(ln, 'SYS / # / OBS TYPES'))
        sys = ln(1);
        nT = str2double(strtrim(ln(2:6)));
        toks = strsplit(strtrim(ln(7:60)));
        toks = toks(~cellfun('isempty', toks));
        if numel(toks) < nT
            pendingSys = sys; pendingTypes = toks; pendingNeed = nT;
        else
            types.(sys) = toks(1:nT);
            order{end+1} = sys; %#ok<AGROW>
        end
    elseif ~isempty(strfind(ln, 'APPROX POSITION XYZ'))
        xyzr = sscanf(ln(1:60), '%f', 3)';
        if numel(xyzr) < 3, xyzr = [xyzr zeros(1, 3-numel(xyzr))]; end
    elseif ~isempty(strfind(ln, 'TIME OF FIRST OBS'))
        v = sscanf(ln(1:60), '%f', 6);
        if numel(v) >= 6, firstObs = v(1:6)'; end
    elseif ~isempty(strfind(ln, 'END OF HEADER'))
        inHeader = false;
    end
end
Obs.XYZr = xyzr;
if isempty(firstObs), firstObs = [2000 1 1 0 0 0]; end
Obs.Date = struct('St', firstObs);
Obs.HeaderSys = order;
if ~isfield(types, 'G')
    fclose(fid);
    error('No GPS (G) observation types in header: %s', obsFullPath);
end
typeG = types.G(:)';
nT = numel(typeG);
Obs.Type = struct('G', {typeG});

% ---------- body ----------
cap = 60000;
epSOD = zeros(1, cap);
prnV  = zeros(1, cap);
datV  = nan(nT, cap);
n = 0;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if isempty(ln) || ln(1) ~= '>', continue; end
    % epoch line: > y mo d h mi sec flag nsat
    e = sscanf(ln(2:end), '%f', 8);
    if numel(e) < 8, continue; end
    sod = e(4)*3600 + e(5)*60 + e(6);
    flag = e(7); nsat = round(e(8));
    if nsat <= 0, continue; end
    if flag > 1
        % special event: skip satellite lines
        for k = 1:nsat, if ~ischar(fgetl(fid)), break; end; end
        continue;
    end
    for k = 1:nsat
        sl = fgetl(fid);
        if ~ischar(sl), break; end
        if isempty(sl), continue; end
        if sl(1) ~= 'G', continue; end  % GPS only
        prn = str2double(strtrim(sl(2:min(3, length(sl)))));
        if isnan(prn), continue; end
        if length(sl) < 3 + nT*16
            sl = [sl repmat(' ', 1, 3 + nT*16 - length(sl))];
        end
        n = n + 1;
        if n > cap
            cap = cap * 2;
            epSOD(1, cap) = 0; prnV(1, cap) = 0; datV(nT, cap) = nan;
        end
        epSOD(n) = sod; prnV(n) = prn;
        for t = 1:nT
            f = sl(4+(t-1)*16 : 3+t*16);
            d = strtrim(f(1:14));
            if ~isempty(d)
                v = str2double(d);
                if ~isnan(v), datV(t, n) = v; end
            end
        end
    end
end
fclose(fid);

if n == 0
    error('No GPS observations parsed from: %s', obsFullPath);
end
% Loss-of-lock written as 0.000: 0.0 code/phase is never physical.
datV = datV(:, 1:n);
for t = 1:nT
    if typeG{t}(1) == 'C' || typeG{t}(1) == 'L'
        datV(t, datV(t, :) == 0) = NaN;
    end
end
Obs.Ep  = struct('G', epSOD(1:n));
Obs.PRN = struct('G', prnV(1:n));
Obs.Data = struct('G', datV);
fprintf('  Parsed GPS obs: %d rows x %d types (%s)\n', n, nT, strjoin(typeG, ' '));
end
