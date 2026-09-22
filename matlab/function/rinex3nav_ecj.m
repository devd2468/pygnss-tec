function Nav = rinex3nav_ecj(navFullPath, sys)
% RINEX3NAV_ECJ  Pure-MATLAB RINEX 3.x navigation reader for Galileo (E),
% BeiDou (C) and QZSS (J). All three use 8-line GPS-like Keplerian records.
%
% Output Nav struct (same convention as rinex3nav_gps):
%   Nav.Eph.(sys) = 34xNrec GPS-layout matrix, Nav.PRN.(sys) = 1xNrec
%
% Column semantics (= GPS layout; only orbit+clock cols drive satpos):
%  E: col32 = BGD(E1,E5a) [assumed 3rd field of line 7], health->0.
%     Week stored as-is (these files use GPS-week numbering).
%  C: col32 = TGD(B1I,B2I) [assumed 3rd field of line 7]. Toe is already on
%     the GPS SOW grid in these files (verified Toe=342000 @ Wed 23:00);
%     BDT-GPS 14 s frame difference therefore needs no correction HERE.
%     If other files store true BDT TOW, add +14 s to Toe (see docs/).
%  J: GPS-identical layout (QZSST ~ GPS time).
%
% satpos_xyz_sbias ignores cols 27-31/33-34 except Ttm; health/accuracy are
% read but never gate the solution (verified), so status-field mappings
% only affect logging, never geometry.

Nav = struct('Com', {{}}, 'Ion', [], 'dTime', [], 'Leap', []);
Nav.Eph = struct();
Nav.PRN = struct();

if ~ismember(sys, {'E', 'C', 'J'})
    error('rinex3nav_ecj:sys', 'System must be E, C or J (got %s).', sys);
end

fid = fopen(navFullPath, 'r');
if fid < 0, error('rinex3nav_ecj:open', 'Cannot open: %s', navFullPath); end

% skip header
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if ~isempty(strfind(ln, 'END OF HEADER')), break; end
end

E = zeros(34, 0); P = zeros(1, 0);
isHdr = @(s) ischar(s) && length(s) >= 9 && ...
    any(s(1) == 'GRECJIS') && ~isempty(regexp(s(1:9), '^[GRECJIS](\d{2}| \d) \d{4} ', 'once'));
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if isempty(strtrim(ln)), continue; end
    if ~isHdr(ln), continue; end
    if ln(1) ~= sys
        continue;  % skip other systems (next loop finds next header)
    end
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
        row = [toc(1)-2000; toc(2); toc(3); toc(4); toc(5); toc(6); ...
               d(1); d(2); d(3); d(4); d(5); d(6); d(7); d(8); d(9); ...
               d(10); d(11); d(12); d(13); d(14); d(15); d(16); d(17); ...
               d(18); d(19); d(20); d(21); d(22); d(23); d(24); d(25); ...
               d(26); d(27); d(28)];
        row(isnan(row)) = 0;
        E(:, end+1) = row; %#ok<AGROW>
        P(:, end+1) = prn; %#ok<AGROW>
    catch
        continue;
    end
end
fclose(fid);

if isempty(P)
    error('rinex3nav_ecj:empty', 'No %s navigation records in: %s', sys, navFullPath);
end
Nav.Eph.(sys) = E;
Nav.PRN.(sys) = P;
fprintf('  Parsed %s nav: %d records, PRNs %s\n', sys, numel(P), mat2str(unique(P)));
end

function v = parseD(s)
if iscell(s), s = s{1}; end
s = strrep(strrep(strtrim(s), 'D', 'E'), 'd', 'e');
if isempty(s)
    v = NaN;
else
    v = str2double(s);
    if isempty(v), v = NaN; end
end
end
