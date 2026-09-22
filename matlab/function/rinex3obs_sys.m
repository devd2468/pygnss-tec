function Obs = rinex3obs_sys(obsFullPath, sys)
% RINEX3OBS_SYS  Pure-MATLAB RINEX 3.x observation reader for ONE system.
% sys: 'E' (Galileo), 'C' (BeiDou), 'J' (QZSS), 'R' (GLONASS, parsed only).
% Same output convention as rinex3obs_gps, keyed by sys:
%   Obs.Type.(sys) 1xM cell, Obs.Data.(sys) MxN, Obs.Ep.(sys) 1xN sod,
%   Obs.PRN.(sys) 1xN, Obs.Date.St, Obs.XYZr, Obs.HeaderSys.
% Long wrapped-looking lines are single physical records (up to ~400 chars
% for 24-type receivers); fgetl handles them directly.

Obs = struct('Com', {{}}, 'XYZr', [0 0 0], 'HeaderSys', {{}});
if ~ismember(sys, {'E', 'C', 'J', 'R'})
    error('rinex3obs_sys:sys', 'System must be E, C, J or R (got %s).', sys);
end

fid = fopen(obsFullPath, 'r');
if fid < 0, error('rinex3obs_sys:open', 'Cannot open: %s', obsFullPath); end

types = struct(); order = {};
xyzr = [0 0 0]; firstObs = [];
pendingTypes = {}; pendingNeed = 0; pendingSys = '';
inHeader = true;
while inHeader
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if length(ln) < 60, ln = [ln repmat(' ', 1, 60-length(ln))]; end
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
        s = ln(1);
        nT = str2double(strtrim(ln(2:6)));
        toks = strsplit(strtrim(ln(7:60)));
        toks = toks(~cellfun('isempty', toks));
        if numel(toks) < nT
            pendingSys = s; pendingTypes = toks; pendingNeed = nT;
        else
            types.(s) = toks(1:nT);
            order{end+1} = s; %#ok<AGROW>
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
Obs.rcvpos = xyzr;  % alias used by TECcalculation_SYS
if isempty(firstObs), firstObs = [2000 1 1 0 0 0]; end
Obs.Date = struct('St', firstObs);
Obs.date = firstObs;  % alias used by TECcalculation_SYS
Obs.HeaderSys = order;
if ~isfield(types, sys)
    fclose(fid);
    error('rinex3obs_sys:nosys', 'No %s observation types in header: %s', sys, obsFullPath);
end
typeS = types.(sys)(:)';
nT = numel(typeS);
Obs.Type = struct(sys, {typeS});

cap = 60000;
epSOD = zeros(1, cap);
prnV  = zeros(1, cap);
datV  = nan(nT, cap);
n = 0;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if isempty(ln) || ln(1) ~= '>', continue; end
    e = sscanf(ln(2:end), '%f', 8);
    if numel(e) < 8, continue; end
    sod = e(4)*3600 + e(5)*60 + e(6);
    flag = e(7); nsat = round(e(8));
    if nsat <= 0, continue; end
    if flag > 1
        for k = 1:nsat, if ~ischar(fgetl(fid)), break; end; end
        continue;
    end
    for k = 1:nsat
        sl = fgetl(fid);
        if ~ischar(sl), break; end
        if isempty(sl), continue; end
        if sl(1) ~= sys, continue; end
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
    error('rinex3obs_sys:empty', 'No %s observations parsed from: %s', sys, obsFullPath);
end
% Loss-of-lock is written as 0.000 by some receivers; 0.0 code/phase is
% never physical, so map exact zeros to NaN (keep Doppler/SNR untouched).
datV = datV(:, 1:n);
for t = 1:nT
    if typeS{t}(1) == 'C' || typeS{t}(1) == 'L'
        datV(t, datV(t, :) == 0) = NaN;
    end
end
Obs.Ep  = struct(sys, epSOD(1:n));
Obs.PRN = struct(sys, prnV(1:n));
Obs.Data = struct(sys, datV(:, 1:n));
fprintf('  Parsed %s obs: %d rows x %d types (%s)\n', sys, n, nT, strjoin(typeS, ' '));
end
