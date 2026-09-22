function G = readGIM(inxFile)
%READGIM  Parse CODE IONEX GIM (*.INX) file.
% Outputs:
%   G.lat      latitudes (1 x nLat)
%   G.lon      longitudes (1 x nLon)
%   G.tec      VTEC maps (nLat x nLon x nMaps) in TECU, NaN = missing
%   G.rms      RMS maps (nLat x nLon x nMaps) in TECU
%   G.epochs   datenum vector of map epochs (UTC)
%   G.scale    exponent scaling (e.g. -1 -> divide by 10)
%
% 9999 = no data (mapped to NaN). Values stored as integers in file,
% divided by 10^(-EXPONENT). EXPONENT typically -1 (values in 0.1 TECU).
    G = struct();
    fid = fopen(inxFile, 'r');
    if fid == -1
        error('readGIM:cannotOpen', 'Cannot open: %s', inxFile);
    end
    try
        text = fread(fid, '*char')';
        fclose(fid);
    catch ME
        fclose(fid);
        error('readGIM:readError', 'Failed to read %s: %s', inxFile, ME.message);
    end

    rawLines = strsplit(text, '\n');
    lines = {};
    for i = 1:numel(rawLines)
        ln = strtrim(rawLines{i});
        if ~isempty(ln)
            lines{end+1} = ln; %#ok<AGROW>
        end
    end

    lat1 = NaN; lat2 = NaN; dlat = NaN;
    lon1 = NaN; lon2 = NaN; dlon = NaN;
    exponent = NaN;
    baseRadius = NaN;
    nMaps = NaN;
    intervalSec = NaN;
    firstEpoch = NaN;
    lastEpoch = NaN;
    for i = 1:numel(lines)
        s = lines{i};
        if contains(s, 'LAT1 / LAT2 / DLAT', 'IgnoreCase', true) || ...
           contains(s, 'LAT1 /LAT2 /DLAT', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 3
                lat1 = nums(1); lat2 = nums(2); dlat = nums(3);
            end
        end
        if contains(s, 'LON1 / LON2 / DLON', 'IgnoreCase', true) || ...
           contains(s, 'LON1 /LON2 /DLON', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 3
                lon1 = nums(1); lon2 = nums(2); dlon = nums(3);
            end
        end
        if contains(s, 'EXPONENT', 'IgnoreCase', true) && ...
           ~contains(s, 'DIFFERENTIAL', 'IgnoreCase', true) && ...
           ~contains(s, 'CODE', 'IgnoreCase', true) && ...
           ~contains(s, 'REFERENCE', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 1
                exponent = nums(1);
            end
        end
        if contains(s, 'BASE RADIUS', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 1, baseRadius = nums(1); end
        end
        if contains(s, 'EPOCH OF FIRST MAP', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 6, firstEpoch = datenum(nums(1:6)); end
        end
        if contains(s, 'EPOCH OF LAST MAP', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 6, lastEpoch = datenum(nums(1:6)); end
        end
        if contains(s, 'INTERVAL', 'IgnoreCase', true) && ...
           ~contains(s, 'ELEVATION', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 1, intervalSec = nums(1); end
        end
        if contains(s, '# OF MAPS IN FILE', 'IgnoreCase', true)
            nums = sscanf(s, '%f')';
            if numel(nums) >= 1, nMaps = nums(1); end
        end
    end

    if isnan(lat1) || isnan(lon1) || isnan(exponent)
        error('readGIM:badHeader', 'Could not parse grid parameters from %s', inxFile);
    end

    nLat = round((lat2 - lat1) / dlat) + 1;
    nLon = round((lon2 - lon1) / dlon) + 1;
    latVec = linspace(lat1, lat2, nLat)';
    lonVec = linspace(lon1, lon2, nLon)';

    scale = 10^exponent;

    headerEnd = find(contains(lines, 'END OF HEADER', 'IgnoreCase', true), 1);
    if isempty(headerEnd)
        error('readGIM:badHeader', 'No END OF HEADER found in %s', inxFile);
    end

    tec = NaN(nLat, nLon, nMaps);
    rms = NaN(nLat, nLon, nMaps);
    epochs = NaN(1, nMaps);

    rowsPerBand = ceil(nLon / 15);

    idx = headerEnd;
    for m = 1:nMaps
        while idx <= numel(lines) && ~contains(lines{idx}, 'START OF TEC MAP', 'IgnoreCase', true)
            idx = idx + 1;
        end
        if idx > numel(lines), break; end
        idx = idx + 1;

        if idx > numel(lines), break; end
        nums = sscanf(lines{idx}, '%f')';
        if numel(nums) >= 6, epochs(m) = datenum(nums(1:6)); end
        idx = idx + 1;

        if idx > numel(lines), break; end
        idx = idx + 1;

        for r = 1:nLat
            if idx > numel(lines), break; end
            idx = idx + 1;

            vals = [];
            for row = 1:rowsPerBand
                if idx > numel(lines), break; end
                rowVals = sscanf(lines{idx}, '%f')';
                vals = [vals, rowVals]; %#ok<AGROW>
                idx = idx + 1;
            end
            if numel(vals) > 0
                nToFill = min(nLon, numel(vals));
                tec(r, 1:nToFill, m) = vals(1:nToFill) * scale;
            end
            tec(r, :, m) = nan2na(tec(r, :, m));
        end
    end

    idx = headerEnd;
    for m = 1:nMaps
        while idx <= numel(lines) && ~contains(lines{idx}, 'START OF RMS MAP', 'IgnoreCase', true)
            idx = idx + 1;
        end
        if idx > numel(lines), break; end
        idx = idx + 1;

        if idx > numel(lines), break; end
        nums = sscanf(lines{idx}, '%f')';
        if numel(nums) >= 6, epochs(m) = max(epochs(m), datenum(nums(1:6))); end
        idx = idx + 1;

        if idx > numel(lines), break; end
        idx = idx + 1;

        for r = 1:nLat
            if idx > numel(lines), break; end
            idx = idx + 1;
            vals = [];
            for row = 1:rowsPerBand
                if idx > numel(lines), break; end
                rowVals = sscanf(lines{idx}, '%f')';
                vals = [vals, rowVals]; %#ok<AGROW>
                idx = idx + 1;
            end
            if numel(vals) > 0
                nToFill = min(nLon, numel(vals));
                rms(r, 1:nToFill, m) = vals(1:nToFill) * scale;
            end
        end
    end

    rms(rms >= 9999) = NaN;

    G.lat = latVec;
    G.lon = lonVec;
    G.tec = tec;
    G.rms = rms;
    G.epochs = epochs;
    G.scale = scale;
    G.lat1 = lat1; G.dlat = dlat;
    G.lon1 = lon1; G.dlon = dlon;
    G.exponent = exponent;
    G.nMaps = nMaps;
    G.nLat = nLat;
    G.nLon = nLon;
end

function x = nan2na(x)
    x(x >= 9999) = NaN;
end