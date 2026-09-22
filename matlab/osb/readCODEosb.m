function osb = readCODEosb(osbFile)
% READCODEOSB  Parse CODE MGEX Bias-SINEX OSB file (satellite code biases).
% Portable pure-MATLAB reader; no binaries needed.
%
% Output: struct with dynamic per-satellite fields, e.g.
%   osb.G01.C1C = -6.5673   (nanoseconds; NaN if absent)
%   osb.meta = struct('file', ..., 'nRec', ...)
% Covers all constellations present (G/R/E/C/J); caller filters by need.
%
% Line format (whitespace-tokenized, robust to blank STATION column):
%   OSB <SVN> <PRN> [STATION] <OBS> <START> <END> <UNIT> <VALUE> <STD>

osb = struct();
osb.meta = struct('file', osbFile, 'nRec', 0);

fid = fopen(osbFile, 'r');
if fid < 0
    error('readCODEosb:open', 'Cannot open OSB file: %s', osbFile);
end
n = 0;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if isempty(ln) || ln(1) == '*' || ln(1) == '+', continue; end
    toks = strsplit(strtrim(ln));
    toks = toks(~cellfun('isempty', toks));
    if numel(toks) < 8 || ~strcmp(toks{1}, 'OSB')
        continue;
    end
    % find unit token; value follows it; obs code precedes date tokens
    ui = find(strcmp(toks, 'ns'), 1);
    if isempty(ui) || ui+1 > numel(toks)
        continue;
    end
    di = find(~cellfun('isempty', regexp(toks, '^\d{4}:', 'once')), 1);
    if isempty(di) || di < 4
        continue;
    end
    prnTok = toks{3};
    if numel(prnTok) < 3
        continue;
    end
    code = toks{di-1};
    bias = str2double(toks{ui+1});
    if isnan(bias)
        continue;
    end
    fld = sprintf('%s%02d', prnTok(1), str2double(prnTok(2:3)));
    if isnan(str2double(prnTok(2:3)))
        continue;
    end
    if ~isfield(osb, fld)
        osb.(fld) = struct();
    end
    osb.(fld).(code) = bias;
    n = n + 1;
end
fclose(fid);
osb.meta.nRec = n;
fprintf('  OSB parsed: %d records from %s\n', n, osbFile);
end
