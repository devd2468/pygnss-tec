function refpos = reffromIPPindex(station_name)
% Read station name from PPP file
% Input: 
%        station_name = station name (string arrays)
% Output:
%        refpos       = reference position (PPP solution)

refpos = [];
% Portable lookup: current folder first, then the repo root (works from any cwd)
cands = {fullfile(pwd, 'PPPindex.txt')};
try
    cands{end+1} = fullfile(repoRoot(), 'PPPindex.txt'); %#ok<AGROW>
catch
end
for k = 1:numel(cands)
    readO = fopen(cands{k}, 'r');
    if readO < 0, continue; end
    try
        text = textscan(readO, '%s', 'Delimiter', '\n');
        fclose(readO);
        text = text{1, 1};
        line1 = find(~cellfun('isempty', strfind(text, station_name)), 1);
        if isempty(line1), continue; end
        all_epoch_line = text{line1, 1};
        all_epoch = regexp(all_epoch_line, '\=*', 'split');
        eval(sprintf('refpos= %s;', all_epoch{2}));
        return;
    catch
        try, fclose(readO); catch, end
    end
end
disp('No PPP reference or wrong station')
