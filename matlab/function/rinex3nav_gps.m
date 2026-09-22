function Nav = rinex3nav_gps(navFullPath)
% RINEX3NAV_GPS  Pure-MATLAB RINEX 3.x navigation reader (GPS only).
% Replaces the ReadNAVrinex304.mexw64 binary for the GPS path.
%
% Output: Nav struct mimicking the MEX output used by readrinex304.m:
%   Nav.Com, Nav.Ion (struct GPSA/GPSB or []), Nav.dTime ([]),
%   Nav.Leap, Nav.Eph (.G = 34xNrec), Nav.PRN (.G = 1xNrec)
%
% eph column layout (must match satpos_xyz_sbias.m):
%  1-6: Toc y(4-digit),mo,d,h,mi,s | 7:clk bias 8:drift 9:driftrate
%  10:IODE 11:Crs 12:Dn 13:M0 14:Cuc 15:e 16:Cus 17:sqrtA 18:Toe
%  19:Cic 20:Omega0 21:Cis 22:i0 23:Crc 24:omega 25:Omegadot 26:IDOT
%  27:CodesL2 28:GPSweek 29:L2Pflag 30:accuracy 31:health 32:TGD
%  33:IODC 34:TransTime

Nav = struct('Com', {{}}, 'Ion', [], 'dTime', [], 'Leap', []);
Nav.Eph = struct('G', zeros(34, 0));
Nav.PRN = struct('G', zeros(1, 0));

fid = fopen(navFullPath, 'r');
if fid < 0, error('Cannot open navigation file: %s', navFullPath); end

gpsa = []; gpsb = [];
% ---------- header ----------
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if length(ln) < 60, ln = [ln repmat(' ', 1, 60-length(ln))]; end
    if ~isempty(strfind(ln, 'IONOSPHERIC CORR'))
        tag = strtrim(ln(1:4));
        v = parseD(ln(5:60));
        if strcmp(tag, 'GPSA'), gpsa = v; elseif strcmp(tag, 'GPSB'), gpsb = v; end
    elseif ~isempty(strfind(ln, 'LEAP SECONDS'))
        Nav.Leap = str2double(strtrim(ln(1:6)));
    elseif ~isempty(strfind(ln, 'END OF HEADER'))
        break;
    end
end
if ~isempty(gpsa) || ~isempty(gpsb)
    Nav.Ion = struct('GPSA', gpsa, 'GPSB', gpsb);
end

% ---------- body ----------
E = zeros(34, 0); P = zeros(1, 0);
isHdr = @(s) ischar(s) && length(s) >= 9 && ...
    any(s(1) == 'GRECJIS') && ~isempty(regexp(s(1:9), '^[GRECJIS](\d{2}| \d) \d{4} ', 'once'));
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if isempty(strtrim(ln)), continue; end
    if ~isHdr(ln), continue; end
    sys = ln(1);
    if sys ~= 'G'
        continue;  % skip non-GPS records (next loop finds next header)
    end
    % GPS record: header + 7 data lines
    lines = cell(1, 8);
    lines{1} = ln;
    ok = true;
    for k = 2:8
        l2 = fgetl(fid);
        if ~ischar(l2), ok = false; break; end
        lines{k} = l2;
    end
    if ~ok, break; end
    try
        prn = str2double(strtrim(lines{1}(2:3)));
        toc = sscanf(lines{1}(4:23), '%f', 6);
        if numel(toc) < 6 || isnan(prn), continue; end
        % data fields are 19-wide: line1 starts [24,43,62], lines 2-8 [5,24,43,62]
        d = zeros(1, 3 + 7*4);
        for k = 1:8
            L = lines{k};
            if length(L) < 80, L = [L repmat(' ', 1, 80-length(L))]; end
            if k == 1
                cols = [24, 43, 62];
                for f = 1:3, d(f) = parseD(L(cols(f):cols(f)+18)); end
            else
                cols = [5, 24, 43, 62];
                for f = 1:4, d(3+(k-2)*4+f) = parseD(L(cols(f):cols(f)+18)); end
            end
        end
        % d layout: [bias drift drate | IODE Crs Dn M0 | Cuc e Cus sqrtA |
        %  Toe Cic Om0 Cis | i0 Crc om OmDot | IDOT CdL2 Week L2P |
        %  acc health TGD IODC | Trans Fit s s]
        row = [toc(1)-2000; toc(2); toc(3); toc(4); toc(5); toc(6); ...
               d(1); d(2); d(3); d(4); d(5); d(6); d(7); d(8); d(9); ...
               d(10); d(11); d(12); d(13); d(14); d(15); d(16); d(17); ...
               d(18); d(19); d(20); d(21); d(22); d(23); d(24); d(25); ...
               d(26); d(27); d(28)];
        row(isnan(row)) = 0;
        % status fields that must be numeric (blank->0 already)
        E(:, end+1) = row; %#ok<AGROW>
        P(:, end+1) = prn; %#ok<AGROW>
    catch
        continue;
    end
end
fclose(fid);

if isempty(P)
    error('No GPS navigation records parsed from: %s', navFullPath);
end
Nav.Eph.G = E;
Nav.PRN.G = P;
fprintf('  Parsed GPS nav: %d records, PRNs %s\n', numel(P), mat2str(unique(P)));
end

function v = parseD(s)
% Parse Fortran D-exponent float(s); blank -> NaN
if iscell(s), s = s{1}; end
s = strrep(strrep(strtrim(s), 'D', 'E'), 'd', 'e');
if isempty(s)
    v = NaN;
else
    v = str2double(s);
    if isempty(v), v = NaN; end
end
end
