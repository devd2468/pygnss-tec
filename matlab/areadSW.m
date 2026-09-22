function [t, v] = areadSW(csvFile, colName)
% AREADSW  Read a space-weather index CSV (timestamp,value) robustly.
% Returns datetime vector t (UTC) and double vector v (NaN for gaps).
    t = []; v = [];
    if ~exist(csvFile, 'file')
        return;
    end
    try
        T = readtable(csvFile, 'TextType', 'string');
    catch
        return;
    end
    if ~any(strcmp(T.Properties.VariableNames, 'timestamp'))
        return;
    end
    if nargin < 2 || isempty(colName)
        % second column by default
        colName = T.Properties.VariableNames{2};
    end
    if ~any(strcmp(T.Properties.VariableNames, colName))
        % readtable sanitizes names (e.g. F10.7 -> F10_7)
        colName = strrep(colName, '.', '_');
    end
    if ~any(strcmp(T.Properties.VariableNames, colName))
        return;
    end
    raws = string(T.timestamp);
    % strip timezone suffix like +00:00
    raws = regexprep(raws, '\+.*$', '');
    try
        t = datetime(raws, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', 'TimeZone', 'UTC');
    catch
        try
            t = datetime(raws, 'TimeZone', 'UTC');
        catch
            return;
        end
    end
    vv = T.(colName);
    if iscell(vv)
        v = str2double(string(vv));
    else
        v = double(vv);
    end
end
