function stationData = loadStationPositions(pppFile)
% Load station positions from PPPindex.txt

    stationData = [];
    if ~exist(pppFile, 'file')
        return;
    end

    fid = fopen(pppFile, 'r');
    text = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    textLines = text{1};

    for i = 1:length(textLines)
        line = textLines{i};
        tokens = regexp(line, '^\s*(\w+)\s*=\s*\[(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\]', 'tokens');
        if ~isempty(tokens)
            stationData(end+1).name = tokens{1}{1};
            stationData(end).xyz = [str2double(tokens{1}{2}), ...
                str2double(tokens{1}{3}), str2double(tokens{1}{4})];
        end
    end
end
