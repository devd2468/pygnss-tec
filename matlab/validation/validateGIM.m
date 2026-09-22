function validateGIM(dayRoots, outDir)
%VALIDATEGIM  Compare pipeline VTEC with CODE IONEX GIM VTEC.
% Output: GIM_validation_stats.csv + GIM_validation.png in outDir
% Note: Pipeline VTEC contains receiver DCB offset, so comparison is for
% shape/timing agreement (correlation) not absolute calibration.
%
% If outDir not provided, uses <dayRoots{1}>/../ArticleFigs
% If dayRoots not provided, uses config dataRoot as single day

    if nargin < 2 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    if nargin < 1 || isempty(dayRoots)
        cfg = config();
        dayRoots = {cfg.dataRoot};
    end

    cfg = config('dataRoot', dayRoots{1});

    % Load station ECEF -> LLA positions
    sp = loadStationPositions(fullfile(repoRoot(), 'PPPindex.txt'));
    if isempty(sp), sp = getDefaultStationData(); end
    stationLLA = containers.Map;
    for i = 1:numel(sp)
        name = char(sp(i).name);
        xyz = sp(i).xyz;
        ll = ecef2lla(xyz);
        stationLLA(name) = struct('lat', ll(1), 'lon', ll(2));
    end

    % Load GIM grid from first file (all same month)
    gimDir = 'D:\IGS DATA\May 2025';
    gimFiles = dir(fullfile(gimDir, 'COD0OPSFIN_*.INX'));
    if isempty(gimFiles), error('validateGIM: No GIM INX files found in %s', gimDir); end
    G0 = readGIM(fullfile(gimDir, gimFiles(1).name));
    nLatG = numel(G0.lat); nLonG = numel(G0.lon);
    latG = G0.lat;
    lonG = G0.lon;

    % Find all stations that have MultiGNSS files across all dayRoots
    stationSet = {};
    for d = 1:numel(dayRoots)
        mfiles = dir(fullfile(dayRoots{d}, 'Results', 'MultiGNSS_*.mat'));
        for k = 1:numel(mfiles)
            m = regexp(mfiles(k).name, '^MultiGNSS_([A-Z0-9]{4})_', 'tokens', 'once');
            if ~isempty(m)
                stn = m{1};
                if ~ismember(stn, stationSet)
                    stationSet{end+1} = stn; %#ok<AGROW>
                end
            end
        end
    end
    stations = unique(stationSet);

    fprintf('validateGIM: Processing %d stations (%s)...\n', numel(stations), strjoin(stations, ', '));

    allStats = struct([]);

    for s = 1:numel(stations)
        stn = stations{s};
        fprintf('  Station %s...\n', stn);

        if ~isKey(stationLLA, stn)
            fprintf('    Station %s not found in PPPindex.txt; skipping.\n', stn);
            continue;
        end

        latS = stationLLA(stn).lat;
        lonS = stationLLA(stn).lon;

        if latS < min(latG) || latS > max(latG) || lonS < min(lonG) || lonS > max(lonG)
            fprintf('    Station %s outside GIM grid (lat=%.2f, lon=%.2f); skipping.\n', ...
                stn, latS, lonS);
            continue;
        end

        [~, ilat] = min(abs(latG - latS));
        [~, ilon] = min(abs(lonG - lonS));
        ilat = max(1, min(nLatG, ilat));
        ilon = max(1, min(nLonG, ilon));

        dayCorr = []; dayBias = []; dayRmse = [];

        for d = 1:numel(dayRoots)
            searchPat = fullfile(dayRoots{d}, 'Results', sprintf('MultiGNSS_%s_*.mat', stn));
            allMFiles = dir(searchPat);
            if isempty(allMFiles), continue; end

            for k = 1:numel(allMFiles)
                fname = allMFiles(k).name;
                m = regexp(fname, '^MultiGNSS_[A-Z0-9]{4}_(\d{4}_\d{2}_\d{2})\.mat$', 'tokens', 'once');
                if isempty(m), continue; end
                ds = m{1};

                try
                    S = load(fullfile(dayRoots{d}, 'Results', fname));
                    V = S.TEC.vertical;
                    pv = fiveMinMed(V);

                    year = str2double(ds(1:4));
                    month = str2double(ds(6:7));
                    day = str2double(ds(9:10));
                    doy = datenum(year, month, 1) + day - 1 - datenum(year, 1, 1) + 1;
                    gimDOY = floor(doy);
                    gimFile = fullfile(gimDir, sprintf('COD0OPSFIN_%d%03d0000_01D_01H_GIM.INX', year, gimDOY));

                    if ~exist(gimFile, 'file')
                        fprintf('      GIM file not found for day %s: %s\n', ds, gimFile);
                        continue;
                    end

                    G = readGIM(gimFile);
                    if isempty(G) || ~isfield(G, 'tec') || isempty(G.tec)
                        fprintf('      Failed to read GIM file: %s\n', gimFile);
                        continue;
                    end

                    [~, ilat] = min(abs(latG - latS));
                    [~, ilon] = min(abs(lonG - lonS));
                    ilat = max(1, min(nLatG-1, ilat));
                    ilon = max(1, min(nLonG-1, ilon));

                    % unit weights: fraction along each axis in [0,1]
                    % (latitude axis may be descending -> signed division is fine)
                    fa = (latS - latG(ilat)) / (latG(ilat+1) - latG(ilat));
                    fb = (lonS - lonG(ilon)) / (lonG(ilon+1) - lonG(ilon));
                    fa = max(0, min(1, fa));
                    fb = max(0, min(1, fb));
                    wA = [1-fa; fa];   % lat weights (rows)
                    wB = [1-fb; fb];   % lon weights (cols)
                    W  = wA * wB';     % 2x2 weight matrix

                    nPipelines = numel(pv);
                    nGimMaps = numel(G.epochs);
                    if nGimMaps < 1, continue; end

                    rDiff  = zeros(1, nPipelines);   % pipeline - GIM
                    gimSer = zeros(1, nPipelines);   % GIM VTEC series
                    pipeSer = zeros(1, nPipelines);  % pipeline VTEC series

                    for k = 1:nPipelines
                        gimHourIdx = floor((k-1) / 12);
                        gimHourIdx = max(1, min(nGimMaps, gimHourIdx));

                        vC = [G.tec(ilat,   ilon,   gimHourIdx), ...
                              G.tec(ilat,   ilon+1, gimHourIdx); ...
                              G.tec(ilat+1, ilon,   gimHourIdx), ...
                              G.tec(ilat+1, ilon+1, gimHourIdx)];

                        gimV = sum(sum(W .* vC));
                        pipeV = pv(k);
                        gimSer(k)  = gimV;
                        pipeSer(k) = pipeV;
                        rDiff(k)   = pipeV - gimV;
                    end

                    % keep samples where BOTH series are finite and sane
                    valid = isfinite(rDiff) & isfinite(gimSer) & isfinite(pipeSer) ...
                            & abs(rDiff) < 30;
                    if sum(valid) < 20
                        fprintf('      Day %s: too few valid points (%d), skipping.\n', ds, sum(valid));
                        continue;
                    end

                    x = pipeSer(valid); y = gimSer(valid);
                    Cmat = corrcoef(x, y);
                    dayCorr(end+1) = Cmat(1, 2);          % pipeline-vs-GIM Pearson r
                    dayBias(end+1) = mean(x - y);
                    dayRmse(end+1) = sqrt(mean(rDiff(valid).^2));

                    fprintf('      Day %s: n=%d, r(pipe,GIM)=%.4f, bias=%.2f TECU, RMSE=%.2f TECU\n', ...
                        ds, sum(valid), dayCorr(end), dayBias(end), dayRmse(end));

                catch ME
                    fprintf('      Day %s error: %s\n', ds, ME.message);
                end
            end
        end

        if ~isempty(dayCorr)
            allStats(end+1).station = stn;
            allStats(end).corr = mean(dayCorr);
            allStats(end).bias = mean(dayBias);
            allStats(end).rmse = mean(dayRmse);
            allStats(end).nDays = numel(dayCorr);
        end
    end

    if ~isempty(allStats) && isfield(allStats(1), 'station')
        fid = fopen(fullfile(outDir, 'GIM_validation_stats.csv'), 'w');
        fprintf(fid, 'station,r_pipeline_GIM,bias_mean_TECU,rmse_TECU,nDays\n');
        for i = 1:numel(allStats)
            if ~isfield(allStats(i), 'station'), continue; end
            fprintf(fid, '%s,%.4f,%.2f,%.2f,%d\n', allStats(i).station, allStats(i).corr, ...
                allStats(i).bias, allStats(i).rmse, allStats(i).nDays);
        end
        fclose(fid);

        fig = figure('Name', 'GIM Validation', 'Position', [100 50 1200 600], 'Visible', 'off');
        nStations = numel(allStats);
        nCols = min(3, nStations);
        nRows = ceil(nStations / nCols);

        for i = 1:nStations
            subplot(nRows, nCols, i);
            hold on; grid on;
            text(0.5, 0.5, sprintf('Station: %s\nBias: %.2f TECU\nRMSE: %.2f TECU\nCorr: %.4f\nDays: %d', ...
                allStats(i).station, allStats(i).bias, allStats(i).rmse, allStats(i).corr, allStats(i).nDays), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9);
            title(allStats(i).station, 'FontSize', 10);
        end

        outPng = fullfile(outDir, 'GIM_validation.png');
        saveas(fig, outPng);
        close(fig);
        fprintf('GIM validation saved: %s\n', outDir);
    else
        fprintf('validateGIM: no data processed.\n');
    end
end

function vm = fiveMinMed(V)
    nb = floor(size(V, 1) / 300);
    vm = nan(nb, 1);
    for b = 1:nb
        blk = V((b-1)*300+1:min(b*300, size(V, 1)), :);
        vm(b) = nanmedian(blk(:));
    end
end