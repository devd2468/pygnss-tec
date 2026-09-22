%% Multi-GNSS TEC Processor - Single Main Entry Point
% ============================================================
% ONE main file covering the full pipeline:
%   1. Multi-constellation RINEX 3.04 intake + TEC: GPS (validated L1/L2),
%      Galileo E1/E5a(E5b), BeiDou B1I/B2I (B1C/B2a fb), QZSS L1/L2(L5);
%      GLONASS detected-only (FDMA/integrator deferred)
%   2. Broadcast nav handling (explicit nav files, no bad auto-download)
%   3. TEC + ROTI + DCB per station (original CSSRG functions)
%   4. Multi-station spatial VTEC mapping (IPP-ready, 350 km shell const)
%   5. ROTI irregularity monitoring + DCB stability + dashboard
%   6. ALL 8 Python plot replications (plot*.m)
% ============================================================
% CONFIG: all machine-specific paths + selections live in config.m
% USAGE:
%   >> MultiGNSS_Main                                     % uses config()
%   >> MultiGNSS_Main('/data/igs/2025-05-01')             % positional dataRoot
%   >> MultiGNSS_Main('dataRoot', X, 'stations', {'BHPL'})% name/value overrides
% ============================================================
function MultiGNSS_Main(varargin)
    inArgs = varargin; %#ok<NASGU>  % capture before clear
    close all; clc
    warning off
    tic_start = tic;

%% ---------- paths ----------
    progRoot = fileparts(mfilename('fullpath'));
    if isempty(progRoot), progRoot = pwd; end
    addpath(progRoot);
    addpath(fullfile(progRoot, 'function'));
    addpath(fullfile(progRoot, 'osb'));    % readCODEosb
    addpath(fullfile(progRoot, 'RINEX'));  % legacy mex/util dir

%% ---------- config (single source: config.m) ----------
    cfgArgs = {};
    if numel(inArgs) >= 1 && ischar(inArgs{1}) && exist(inArgs{1}, 'dir')
        cfgArgs = [{'dataRoot', inArgs{1}}];           % positional data dir
    elseif ~isempty(inArgs)
        cfgArgs = inArgs;                              % name/value overrides
    end
    cfg = config(cfgArgs{:});
    DATA_ROOT = cfg.dataRoot;
    MGEX_ROOT = cfg.mgexRoot;
    OUTPUT_ROOT = cfg.outRoot;
    CONSTELLATIONS = cfg.constellations;
    if ~exist(OUTPUT_ROOT, 'dir'), mkdir(OUTPUT_ROOT); end
    if ~isempty(cfg.stations)
        fprintf('Station filter: %s\n', strjoin(cfg.stations, ', '));
    end

    fprintf('\n====================================================\n');
    fprintf('  MULTI-GNSS TEC PROCESSOR - COMPLETE PIPELINE\n');
    fprintf('  Constellations: %s\n', strjoin(CONSTELLATIONS, ', '));
    fprintf('====================================================\n');
    fprintf('Program root : %s\n', progRoot);
    fprintf('Data dir     : %s\n', DATA_ROOT);
    fprintf('MGEX dir     : %s\n', MGEX_ROOT);
    fprintf('Output dir   : %s\n\n', OUTPUT_ROOT);

%% ---------- binaries (mex/curl/gzip) must be visible from DATA dir ----------
    ensureBinaries(progRoot, DATA_ROOT);

%% ---------- MGEX inventory (staged for precise-orbit use) ----------
    if exist(MGEX_ROOT, 'dir')
        sp3 = dir(fullfile(MGEX_ROOT, '*_ORB.SP3'));
        clk = dir(fullfile(MGEX_ROOT, '*_CLK.CLK'));
        erp = dir(fullfile(MGEX_ROOT, '*_ERP.ERP'));
        osb = dir(fullfile(MGEX_ROOT, '*_OSB.BIA'));
        fprintf('MGEX products: SP3=%d CLK=%d ERP=%d OSB=%d\n\n', ...
            numel(sp3), numel(clk), numel(erp), numel(osb));
    else
        fprintf('MGEX dir not found, continuing with broadcast ephemeris.\n\n');
    end

%% ---------- list OBSERVATION files only (never navigation) ----------
    obsFiles = listObsFiles(DATA_ROOT);
    if ~isempty(cfg.stations)
        keep = false(numel(obsFiles), 1);
        for i = 1:numel(obsFiles)
            keep(i) = any(strncmp(obsFiles(i).name, cfg.stations, 4));
        end
        obsFiles = obsFiles(keep);
    end
    if isempty(obsFiles)
        error('No RINEX observation files (*_MO.rnx / *.??o) in: %s', DATA_ROOT);
    end
    fprintf('Found %d observation file(s) (nav files excluded):\n', numel(obsFiles));
    for i = 1:numel(obsFiles)
        fprintf('  %s (%.1f MB)\n', obsFiles(i).name, obsFiles(i).bytes/1e6);
    end
    fprintf('\n');

%% ---------- process each station ----------
    allResults = {};
    nOK = 0;
    for i = 1:numel(obsFiles)
        r_o_name = obsFiles(i).name;
        station_name = r_o_name(1:4);
        fprintf('\n====================================================\n');
        fprintf('  STATION %s : %s\n', station_name, r_o_name);
        fprintf('====================================================\n');
        try
            navName = findNavFile(DATA_ROOT, r_o_name);
            if ~isempty(navName)
                fprintf('  Navigation file: %s\n', navName);
            else
                % No local nav + legacy FTP auto-download is unreachable:
                % skip fast instead of hanging on curl timeouts.
                fprintf('  SKIP station %s: no local nav file (no download in offline mode).\n', station_name);
                continue;
            end

            checkfileRN(r_o_name, [DATA_ROOT filesep]);
            fprintf('  Reading RINEX...\n');
            [obs, nav, doy, Year] = readrinex304(r_o_name, navName, [DATA_ROOT filesep]);
            fprintf('  Date %d-%02d-%02d | epochs %d | station %s\n', ...
                obs.date(1), obs.date(2), obs.date(3), numel(obs.epoch), obs.station);
            if isfield(obs, 'headerSys') && ~isempty(obs.headerSys)
                fprintf('  Header constellations: %s\n', strjoin(obs.headerSys, ','));
            end
            if isfield(obs, 'sys')
                extra = fieldnames(obs.sys);
                if ~isempty(extra)
                    fprintf('  Extra constellation data present: %s (staged)\n', strjoin(extra', ','));
                end
            end

            DCB_path = fullfile(DATA_ROOT, 'DCB');
            if ~exist(DCB_path, 'dir'), mkdir(DCB_path); end
            ensureBinaries(progRoot, DCB_path);
            % Satellite DCB priority: local MGEX OSB -> monthly DCB -> zeros.
            % (Legacy CODE FTP is unreachable; never attempt downloads.)
            [satb.P1C1, satb.P1P2, dcbSrc] = getSatDCB(obs, MGEX_ROOT, DCB_path, progRoot);
            fprintf('  Satellite DCB source: %s\n', dcbSrc);

            fprintf('  Calculating TEC (GPS L1/L2 validated path)...\n');
            TECcalculationRINEX304_OEM7(obs, nav, satb, [OUTPUT_ROOT filesep]);

            year = num2str(obs.date(1)); month = num2str(obs.date(2), '%.2d'); date = num2str(obs.date(3), '%.2d');
            resFile = fullfile(OUTPUT_ROOT, sprintf('TEC_%s_%s_%s_%s.mat', obs.station, year, month, date));
            S = load(resFile);
            TEC  = S.(sprintf('TEC_%s_%s_%s',  year, month, date));
            DCB  = S.(sprintf('DCB_%s_%s_%s',  year, month, date));
            ROTI = S.(sprintf('ROTI_%s_%s_%s', year, month, date));
            prm  = S.(sprintf('prm_%s_%s_%s',  year, month, date));
            % Load S4 — saved as S4_<station>_<date> by TEC engine
            s4varname = sprintf('S4_%s_%s_%s_%s', obs.station, year, month, date);
            if isfield(S, s4varname)
                S4 = S.(s4varname);
            else
                S4 = nan(size(ROTI));
                fprintf('  S4 variable not found in %s (all-NaN).\n', resFile);
            end

            outName = fullfile(OUTPUT_ROOT, sprintf('MultiGNSS_%s_%s_%s_%s.mat', obs.station, year, month, date));
            save(outName, 'TEC', 'ROTI', 'DCB', 'prm', 'S4', 'station_name', 'obs', 'nav');
            fprintf('  Saved: %s\n', outName);

            allResults{end+1} = struct('station', obs.station, 'TEC', TEC, ...
                'ROTI', ROTI, 'DCB', DCB, 'prm', prm, 'S4', S4); %#ok<AGROW>
            % ---- Galileo / BeiDou / QZSS (same-station, same-day) ----
            sysRes = struct();
            for sysCell = {'E', 'C', 'J'}
                sys = sysCell{1};
                if ~ismember(sys, CONSTELLATIONS), continue; end
                try
                    r = processSystemSYS(sys, r_o_name, navName, DATA_ROOT, ...
                        OUTPUT_ROOT, MGEX_ROOT, obs);
                    if ~isempty(r)
                        sysRes.(sys) = r;
                    end
                catch ME0
                    fprintf('  %s: skipped (%s)\n', sys, ME0.message);
                end
            end
            if ismember('R', CONSTELLATIONS)
                fprintf('  R (GLONASS): detected-only; FDMA/integrator out of scope.\n');
            end
            if ~isempty(fieldnames(sysRes))
                allResults{end}.sys = sysRes;
            end
            nOK = nOK + 1;
            fprintf('  PASS station %s\n', obs.station);
        catch ME
            fprintf('  FAIL station %s: %s\n', station_name, ME.message);
            for k = 1:min(3, numel(ME.stack))
                fprintf('    at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
    end
    fprintf('\nProcessed %d of %d station(s).\n', nOK, numel(obsFiles));
    if nOK == 0
        error('No stations processed successfully. See errors above.');
    end

%% ---------- spatial mapping ----------
    fprintf('\n====================================================\n  SPATIAL MAPPING\n====================================================\n');
    if numel(allResults) >= 2
        try
            MultiStationSpatialMap(allResults, OUTPUT_ROOT);
        catch ME
            fprintf('  Spatial mapping failed (non-fatal): %s\n', ME.message);
        end
    else
        fprintf('  Need >=2 stations (have %d); skipping.\n', numel(allResults));
    end

%% ---------- ROTI monitor / DCB / ML classifier / dashboard ----------
    fprintf('\n====================================================\n  ROTI MONITOR + DCB + ML + DASHBOARD\n====================================================\n');
    try, MultiGNSS_ROTI_Monitor(allResults, OUTPUT_ROOT);
    catch ME, fprintf('  ROTI monitor failed (non-fatal): %s\n', ME.message); end
    try, DCB_Stability_Analysis(allResults, OUTPUT_ROOT);
    catch ME, fprintf('  DCB analysis failed (non-fatal): %s\n', ME.message); end
    try, MultiGNSS_MLClassifier(allResults, OUTPUT_ROOT);
    catch ME, fprintf('  ML classifier failed (non-fatal): %s\n', ME.message); end
    try
        MultiGNSS_Dashboard(allResults, OUTPUT_ROOT);
        fprintf('  Dashboard done.\n');
    catch ME
        fprintf('  Dashboard failed: %s\n', ME.message);
    end

%% ---------- 8 python-plot replications ----------
    fprintf('\n====================================================\n  PLOTS (8)\n====================================================\n');
    runPlot('1/8 plotMagneticEquator', OUTPUT_ROOT, ...
        @() plotMagneticEquatorDemo(OUTPUT_ROOT));
    runPlot('2/8 plotMappingRange', OUTPUT_ROOT, ...
        @() plotMappingRangeDemo(OUTPUT_ROOT));
    runPlot('3/8 plotNortheastReceivers', OUTPUT_ROOT, ...
        @() plotNortheastReceivers(OUTPUT_ROOT));
    runPlot('4/8 plotROTIdtrend', OUTPUT_ROOT, ...
        @() plotROTIdtrend(OUTPUT_ROOT));
    runPlot('5/8 plotROTIinDifferentsSites', OUTPUT_ROOT, ...
        @() plotROTIinDifferentsSites(OUTPUT_ROOT));
    runPlot('6/8 plotROTImaps', OUTPUT_ROOT, ...
        @() plotROTImaps(OUTPUT_ROOT));
    runPlot('7/8 plotScatterAndInterpolation', OUTPUT_ROOT, ...
        @() plotScatterAndInterpolation(OUTPUT_ROOT));
    runPlot('8/8 plotMapWithIPP', OUTPUT_ROOT, ...
        @() plotMapWithIPP(OUTPUT_ROOT));
    % S4 scintillation figure (article supplement)
    runPlot('S4 aplotS4', OUTPUT_ROOT, ...
        @() aplotS4({DATA_ROOT}, cfg.stations, OUTPUT_ROOT, ...
                    'titleStr', 'S4 Amplitude Scintillation', 'saveFig', true));

    %% ================================================================
    %% NOVEL RESEARCH ANALYSIS SUITE (Q1 GPS Solutions contribution)
    %% Modules: N1 EPBPropagation | N2 ScintCoupling | N3 TIDSpectral
    %%          N4 ConstelDiversity | N5 SolarMax
    %% All outputs → OUTPUT_ROOT; reproducible (rng=42, no toolboxes)
    %% ================================================================
    if isfield(cfg, 'runNovelAnalysis') && cfg.runNovelAnalysis
        fprintf('\n====================================================\n');
        fprintf('  NOVEL ANALYSIS SUITE — Q1 research modules\n');
        fprintf('====================================================\n');
        try
            % Build dayRoots cell: one entry per unique date found in DATA_ROOT
            dayRootsCells = buildDayRoots(DATA_ROOT);
            if isempty(dayRootsCells)
                dayRootsCells = {DATA_ROOT};   % fallback: single root
            end
            % Pass F10.7 / Kp from cfg if available
            novOpts = {};
            if isfield(cfg,'F107') && ~isempty(cfg.F107)
                novOpts = [novOpts, 'F107', cfg.F107];
            end
            if isfield(cfg,'Kp') && ~isempty(cfg.Kp)
                novOpts = [novOpts, 'Kp', cfg.Kp];
            end
            novelResults = runNovelAnalysis(dayRootsCells, cfg.stations, OUTPUT_ROOT, ...
                'IST', 5.5, 'nBoot', 500, 'rngSeed', 42, novOpts{:});
            fprintf('  Novel analysis: %d errors\n', numel(novelResults.errors));
            fprintf('  Report: %s\n', novelResults.reportPath);
        catch ME_novel
            fprintf('  Novel analysis FAILED: %s\n', ME_novel.message);
        end
    else
        fprintf('\n[Info] Novel analysis skipped (cfg.runNovelAnalysis=false).\n');
        fprintf('  To enable: MultiGNSS_Main(''runNovelAnalysis'',true)\n');
    end

    fprintf('\n====================================================\n  COMPLETE: %d station(s), %.1f sec\n  Output: %s\n====================================================\n', ...
        nOK, toc(tic_start), OUTPUT_ROOT);
    warning on;
end

%% ================= local helpers (unique names only) =================
function runPlot(tag, ~, fn)
    try
        fprintf('  %s\n', tag);
        fn();
        fprintf('    PASS\n');
    catch ME
        fprintf('    FAIL: %s\n', ME.message);
    end
end

function dayRoots = buildDayRoots(dataRoot)
% BUILDDAYROOTS  Discover per-day sub-directories or return dataRoot as-is.
% Looks for sub-directories matching YYYY-MM-DD or YYYYDDD patterns.
    dayRoots = {};
    % Try ISO date subdirs first
    subDirs = dir(dataRoot);
    for k = 1:numel(subDirs)
        if ~subDirs(k).isdir, continue; end
        nm = subDirs(k).name;
        if regexp(nm, '^\d{4}-\d{2}-\d{2}$') || regexp(nm, '^\d{4}\d{3}$')
            dayRoots{end+1} = fullfile(dataRoot, nm); %#ok<AGROW>
        end
    end
    if isempty(dayRoots)
        % No date subdirs found — treat dataRoot itself as single day
        dayRoots = {dataRoot};
    else
        dayRoots = sort(dayRoots);  % chronological order
    end
end

function plotMagneticEquatorDemo(OUTPUT_ROOT)
    figure('Name', 'Magnetic Equator / EIA Crests', 'Visible', 'off');
    plotMagneticEquator(gca, '2021', [-30, 30], [-100, 100]);
    title('Equatorial Ionization Anomaly Model (2021)');
    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotMagneticEquator_Matlab.png'));
    close;
end

function plotMappingRangeDemo(OUTPUT_ROOT)
    figure('Name', 'Mapping Range', 'Visible', 'off');
    plotMappingRange(gca, 100.5, 13.75, 500);
    title('Station Observation Range (500 km)');
    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotMappingRange_Matlab.png'));
    close;
end

function ensureBinaries(progRoot, destDir)
% Legacy helper, now a no-op: pure-MATLAB parsers are default (no MEX) and
% downloads use built-in websave/gunzip (no curl/gzip exes). Kept as a stub
% so older scripts calling it keep working on any OS.
    %#ok<INUSD>
end

function obsFiles = listObsFiles(DATA_ROOT)
% Observation files ONLY: *_MO.rnx or *.??o ; exclude *MN/*EN/*GN/*RN/*SN/*IN
    c1 = dir(fullfile(DATA_ROOT, '*_MO.rnx'));
    c2 = dir(fullfile(DATA_ROOT, '*.??o'));
    c3 = dir(fullfile(DATA_ROOT, '*_O.rnx'));
    all = [c1; c2; c3];
    keep = true(numel(all), 1);
    for k = 1:numel(all)
        nm = all(k).name;
        if ~isempty(regexpi(nm, '_[MEGRSI]N\.rnx$', 'once')) || ...
           ~isempty(regexpi(nm, '\.[0-9]{2}n$', 'once'))
            keep(k) = false;  % navigation file, never feed to obs reader
        end
    end
    obsFiles = all(keep);
    % de-duplicate + sort
    [~, u] = unique({obsFiles.name}, 'stable');
    obsFiles = obsFiles(u);
    [~, o] = sort({obsFiles.name});
    obsFiles = obsFiles(o);
end

function navName = findNavFile(DATA_ROOT, obsName)
% Derive navigation filename from observation filename (basename only,
% because legacy reader cd()s into the data dir).
    navName = '';
    cands = {};
    % IGS long names: ..._30S_MO.rnx -> ..._MN.rnx
    t = regexprep(obsName, '_\d+S_MO\.rnx$', '_MN.rnx', 'ignorecase');
    if ~strcmp(t, obsName), cands{end+1} = t; end
    t = regexprep(obsName, '_MO\.rnx$', '_MN.rnx', 'ignorecase');
    if ~strcmp(t, obsName), cands{end+1} = t; end
    % Classic: .23o -> .23n (any .??o -> .??n)
    t = regexprep(obsName, '\.([0-9A-Za-z]{2})o$', '.$1n');
    if ~strcmp(t, obsName), cands{end+1} = t; end
    % Prefix fallback: first 9 chars (IGS) or 4 chars, any *MN.rnx
    for L = [9, 4]
        if numel(obsName) >= L
            d = dir(fullfile(DATA_ROOT, [obsName(1:L) '*MN.rnx']));
            for k = 1:numel(d), cands{end+1} = d(k).name; end
        end
    end
    for k = 1:numel(cands)
        if exist(fullfile(DATA_ROOT, cands{k}), 'file')
            navName = cands{k};
            return;
        end
    end
end

function [P1C1, P1P2, src] = getSatDCB(obs, MGEX_ROOT, DCB_path, progRoot)
% Satellite DCB priority: local MGEX OSB -> monthly DCB files -> zeros.
% OSB path uses the same frequency pair as the TEC engine (pickFreqPair)
% with P1Cx = OSB_C1W - OSB_Cx (seconds); missing sats stay 0 and counted.
    P1C1 = zeros(32, 1); P1P2 = zeros(32, 1);
    src = 'zeros (no DCB)';
    try
        [tC1, tP2] = pickFreqPair(obs);
        dn = datenum(obs.date(1:3));
        doy = round(dn - datenum(obs.date(1), 1, 1)) + 1;
        key = sprintf('%04d%03d', obs.date(1), doy);
        o = osbForDate(MGEX_ROOT, obs.date(1:3));
        if ~isempty(o)
            nHit = 0;
            for prn = 1:32
                f = sprintf('G%02d', prn);
                if isfield(o, f)
                    s = o.(f);
                    if isfield(s, 'C1W') && isfield(s, tC1) && isfield(s, tP2)
                        P1C1(prn) = (s.C1W - s.(tC1)) * 1e-9;
                        P1P2(prn) = (s.C1W - s.(tP2)) * 1e-9;
                        nHit = nHit + 1;
                    end
                end
            end
            src = sprintf('MGEX OSB day %s pair %s/%s (%d/32 sats)', ...
                key, tC1, tP2, nHit);
            return;
        end
    catch ME
        fprintf('  OSB DCB note: %s\n', ME.message);
    end
    try
        yy = num2str(obs.date(1)); yy = yy(3:4);
        mm = num2str(obs.date(2), '%.2d');
        if exist(fullfile(DCB_path, ['P1C1' yy mm '.DCB']), 'file') && ...
           exist(fullfile(DCB_path, ['P1P2' yy mm '.DCB']), 'file')
            [P1C1, P1P2] = dlsat(obs, progRoot, [DCB_path filesep]);
            src = 'local monthly DCB';
        end
    catch
    end
end

function r = processSystemSYS(sys, r_o_name, navName, DATA_ROOT, OUTPUT_ROOT, MGEX_ROOT, obsGPS)
% One non-GPS constellation for one station-day. Returns struct with
% TEC/ROTI/DCB/prm or [] if unavailable. Never throws (caller guards too).
    r = [];
    if ~ismember(sys, {'E', 'C', 'J'})
        return;
    end
    MO = fullfile(DATA_ROOT, r_o_name);
    MN = fullfile(DATA_ROOT, navName);
    O = rinex3obs_sys(MO, sys);
    N = rinex3nav_ecj(MN, sys);
    % best-validity pair across receiver signal flavors (shared picker)
    [PT, F, best] = pickSysPair(O.Type.(sys), O.Data.(sys)', sys);
    if best == 0 || isempty(PT.c1)
        fprintf('  %s: no complete dual-frequency pair in file.\n', sys);
        return;
    end
    fprintf('  %s pair: %s/%s + %s/%s (k=%.4f, n=%d)\n', ...
        sys, PT.c1, PT.l1, PT.c2, PT.l2, F.k, best);
    % remap to engine struct (NxM data, Nx1 index/epoch)
    obsS = struct('type', {O.Type.(sys)}, 'data', O.Data.(sys)', ...
        'index', O.PRN.(sys)', 'epoch', O.Ep.(sys)', ...
        'rcvpos', O.rcvpos, 'station', obsGPS.station, 'date', obsGPS.date);
    navS = struct('eph', N.Eph.(sys)', 'index', N.PRN.(sys)');
    % satellite DCB from cached OSB (seconds, code2-code1)
    o = osbForDate(MGEX_ROOT, obsGPS.date(1:3));
    sats = unique(obsS.index); sats = sats(sats >= 1 & sats <= 64);
    if ~isempty(o)
        [satBias_s, nHit] = osbPairDCB(o, sys, sats, PT.c1, PT.c2);
        fprintf('  %s sat DCB: OSB %s/%s (%d sats)\n', sys, PT.c1, PT.c2, nHit);
    else
        satBias_s = zeros(64, 1);
        fprintf('  %s sat DCB: zeros (no OSB)\n', sys);
    end
    tag = sprintf('%s_%s', obsGPS.station, sys);
    TECcalculation_SYS(obsS, navS, satBias_s, [OUTPUT_ROOT filesep], tag, F, PT);
    year = num2str(obsGPS.date(1)); month = num2str(obsGPS.date(2), '%.2d');
    date = num2str(obsGPS.date(3), '%.2d');
    resFile = fullfile(OUTPUT_ROOT, sprintf('TEC_%s_%s_%s_%s.mat', tag, year, month, date));
    S = load(resFile);
    T2 = S.(sprintf('TEC_%s_%s_%s_%s', tag, year, month, date));
    D2 = S.(sprintf('DCB_%s_%s_%s_%s', tag, year, month, date));
    R2 = S.(sprintf('ROTI_%s_%s_%s_%s', tag, year, month, date));
    P2 = S.(sprintf('prm_%s_%s_%s_%s', tag, year, month, date));
    s4varname_sys = sprintf('S4_%s_%s_%s_%s', tag, year, month, date);
    if isfield(S, s4varname_sys)
        S4_2 = S.(s4varname_sys);
    else
        S4_2 = nan(size(R2));
    end
    outName = fullfile(OUTPUT_ROOT, sprintf('MultiGNSS_%s_%s_%s_%s_%s.mat', ...
        obsGPS.station, sys, year, month, date));
    TEC = T2; ROTI = R2; DCB = D2; prm = P2; S4 = S4_2; %#ok<NASGU>
    save(outName, 'TEC', 'ROTI', 'DCB', 'prm', 'S4');
    fprintf('  Saved: %s\n', outName);
    r = struct('TEC', T2, 'ROTI', R2, 'DCB', D2, 'prm', P2, 'S4', S4_2, 'pair', PT);
end

function MultiStationSpatialMap(allResults, OUTPUT_ROOT)
    numStations = numel(allResults);
    stationNames = cellfun(@(s) s.station, allResults, 'UniformOutput', false);
    % Station positions: standalone loader returns struct array
    sp = loadStationPositions(fullfile(repoRoot(), 'PPPindex.txt'));
    if isempty(sp), sp = getDefaultStationData(); end
    % Map requested stations to coordinates (ECEF -> LLA)
    xyz = nan(numStations, 3);
    for s = 1:numStations
        ix = find(strcmp({sp.name}, stationNames{s}), 1);
        if isempty(ix)
            ix = find(strncmp({sp.name}, stationNames{s}, 4), 1);
        end
        if ~isempty(ix), xyz(s, :) = sp(ix).xyz; end
    end
    ok = all(~isnan(xyz), 2);
    if sum(ok) < 2
        fprintf('  Spatial map skipped: <2 known station positions.\n');
        return;
    end
    xyz = xyz(ok, :); stationNames = stationNames(ok);
    allResults = allResults(ok);
    numStations = sum(ok);
    lla = zeros(numStations, 3);
    for s = 1:numStations, lla(s, :) = ecef2lla(xyz(s, :)); end
    lats = lla(:, 1); lons = lla(:, 2);

    lat_min = min(lats)-1; lat_max = max(lats)+1;
    lon_min = min(lons)-1; lon_max = max(lons)+1;
    [lonGrid, latGrid] = meshgrid(lon_min:0.1:lon_max, lat_min:0.1:lat_max);

    % Midday VTEC per station, interpolated over grid
    vals = nan(numStations, 1);
    for s = 1:numStations
        try
            V = allResults{s}.TEC.vertical;
            idx = min(size(V, 1), round(12*3600)+1);
            vals(s) = nanmean(V(idx, :));
        catch
        end
    end
    F = griddata(lons, lats, vals, lonGrid, latGrid, 'natural');
    figure('Name', 'VTEC spatial map', 'Visible', 'off');
    pcolor(lonGrid, latGrid, F); shading interp; hold on;
    plot(lons, lats, 'k*', 'MarkerSize', 12);
    for s = 1:numStations
        text(lons(s)+0.15, lats(s)+0.15, stationNames{s}, 'FontSize', 8);
    end
    colorbar; caxis([0 80]);
    xlabel('Longitude'); ylabel('Latitude');
    title('VTEC Map at 12:00 UTC (IDW/natural interpolation)');
    grid on;
    outPng = fullfile(OUTPUT_ROOT, 'MultiGNSS_VTEC_SpatialMap.png');
    try
        saveas(gcf, outPng);
    catch
        try  % headless fallback: explicit painters print
            print(gcf, outPng, '-dpng', '-r150', '-painters');
        catch ME2
            fprintf('  Spatial map export failed (non-fatal): %s\n', ME2.message);
            try, close; catch, end
            return;
        end
    end
    close;
    fprintf('  Spatial map saved.\n');
end

function MultiGNSS_ROTI_Monitor(allResults, OUTPUT_ROOT)
    fprintf('  ROTI monitor: thresholds 0.2 / 0.5 / 1.0 TECU/min\n');
    fprintf('  S4  monitor: thresholds 0.15 (weak) / 0.30 (moderate) / 0.50 (strong)\n');
    evAll = {};
    for s = 1:numel(allResults)
        try
            R = allResults{s}.ROTI;
            med = nanmedian(R, 2);           % 86400 x 1, per-second rows
            tmin = (0:numel(med)-1)/60;      % minutes
            isEv = med >= 0.5;
            % group contiguous ROTI events
            d = diff([0; isEv; 0]);
            st_idx = find(d == 1); en_idx = find(d == -1) - 1;
            nev = 0;
            for k = 1:numel(st_idx)
                dur = tmin(en_idx(k)) - tmin(st_idx(k));
                if dur >= 15, nev = nev + 1; end
            end

            % S4 statistics
            maxS4  = NaN;
            s4nMod = NaN;
            s4nStr = NaN;
            if isfield(allResults{s}, 'S4') && ~isempty(allResults{s}.S4)
                S4mat  = allResults{s}.S4;
                s4med  = nanmedian(S4mat, 2);  % per-second median across PRNs
                maxS4  = nanmax(s4med);
                s4nMod = sum(s4med >= 0.3, 'omitnan');   % epochs >= moderate
                s4nStr = sum(s4med >= 0.5, 'omitnan');   % epochs >= strong
            end

            evAll{end+1} = struct('station', allResults{s}.station, ...
                'maxROTI', nanmax(med), 'nEvents', nev, ...
                'maxS4', maxS4, 's4nModerate', s4nMod, 's4nStrong', s4nStr); %#ok<AGROW>
            fprintf('    %s: maxROTI=%.3f, events(>=0.5,>=15min)=%d | maxS4=%.3f, epochs>=0.3=%d, >=0.5=%d\n', ...
                allResults{s}.station, nanmax(med), nev, maxS4, s4nMod, s4nStr);
        catch ME
            fprintf('    %s ROTI/S4 failed: %s\n', allResults{s}.station, ME.message);
        end
    end
    try, save(fullfile(OUTPUT_ROOT, 'ROTI_events.mat'), 'evAll'); catch, end
end

function DCB_Stability_Analysis(allResults, OUTPUT_ROOT)
    st = cellfun(@(s) s.station, allResults, 'UniformOutput', false);
    rcv = nan(numel(allResults), 1);
    for s = 1:numel(allResults)
        try, rcv(s) = allResults{s}.DCB.rcv; catch, end
    end
    fprintf('  Receiver DCB (TECU): ');
    for s = 1:numel(st), fprintf('%s=%.2f ', st{s}, rcv(s)); end
    fprintf('\n');
    try
        figure('Name', 'DCB stability', 'Visible', 'off');
        bar(rcv); set(gca, 'XTickLabel', st, 'XTickLabelRotation', 45);
        ylabel('Receiver DCB (TECU)'); title('Receiver DCB per station');
        grid on;
        saveas(gcf, fullfile(OUTPUT_ROOT, 'DCB_Stability_Analysis.png'));
        close;
        save(fullfile(OUTPUT_ROOT, 'DCB_Stability_Analysis.mat'), 'st', 'rcv');
    catch ME
        fprintf('  DCB plot failed: %s\n', ME.message);
    end
end

function MultiGNSS_MLClassifier(allResults, OUTPUT_ROOT)
% Phase 3: threshold + Random-Forest disturbance classification (guarded).
% Classes: 0 Quiet (<0.5), 1 Mild EPB (0.5-1.0), 2 Severe (>=1.0 max ROTI).
    try
        F = []; L = [];
        for s = 1:numel(allResults)
            try
                R = allResults{s}.ROTI;
                V = allResults{s}.TEC.vertical;
                E = allResults{s}.prm.elevation;
                nW = floor(86400/300);
                for w = 1:nW
                    i0 = (w-1)*300+1; i1 = min(w*300, 86400);
                    ws = V(i0:i1, :);
                    wr = R(i0:min(size(R,1),i1), :);  % ROTI per-second rows
                    if sum(~isnan(ws(:))) < 10, continue; end
                    f1 = nanmean(wr(:)); f2 = nanmax(wr(:));
                    % 30 s-lag gradient (matrix is per-second, obs every 30 s)
                    w30 = ws(1:30:end, :);
                    d30 = abs(diff(w30, 1, 1)); d30 = d30(~isnan(d30));
                    if isempty(d30), f4 = 0; else, f4 = nanmean(d30); end
                    f5 = nanstd(ws(:))^2;
                    if isnan(f5), f5 = 0; end
                    we = E(i0:i1, :); m = we > 30;
                    f6 = nansum(ws(:).*m(:)) / max(sum(m(:)), 1);
                    f7 = sum(~all(isnan(ws), 1));
                    F(end+1, :) = [f1, f2, f4, f5, f6, f7]; %#ok<AGROW>
                    if f2 >= 1.0, L(end+1) = 2; elseif f2 >= 0.5, L(end+1) = 1; else, L(end+1) = 0; end %#ok<AGROW>
                end
            catch
            end
        end
        ok = all(~isnan(F), 2);
        F = F(ok, :); L = L(ok); L = L(:);  % force column: avoids N-by-N == expansion
        fprintf('  ML features: %d windows (Q=%d M=%d S=%d)\n', ...
            numel(L), sum(L==0), sum(L==1), sum(L==2));
        if size(F, 1) < 20
            fprintf('  ML skipped: insufficient windows.\n');
            return;
        end
        try
            model = TreeBagger(100, F, L, 'Method', 'classification', ...
                'OOBPrediction', 'on', 'MinLeafSize', 5);
            pred = str2double(predict(model, F)); pred = pred(:);
            acc = mean(pred == L);
            fprintf('  Random Forest accuracy: %.2f%%\n', acc*100);
            try
                imp = predictorImportance(model);
            catch
                try
                    imp = model.OOBPermutedPredictorDeltaError;
                catch
                    imp = ones(1, size(F, 2));
                end
            end
            names = {'MeanROTI','MaxROTI','dSTECdt','PhaseVar','ElevSTEC','SatCount'};
            [si, ix] = sort(imp, 'descend');
            for k = 1:numel(ix), fprintf('    %s: %.4f\n', names{ix(k)}, si(k)); end
            figure('Name', 'ML importance', 'Visible', 'off');
            bar(si); set(gca, 'XTickLabel', names(ix), 'XTickLabelRotation', 45);
            ylabel('Importance'); title('RF Feature Importance'); grid on;
            saveas(gcf, fullfile(OUTPUT_ROOT, 'Phase3_FeatureImportance.png'));
            close;
            save(fullfile(OUTPUT_ROOT, 'Phase3_MLModel.mat'), 'acc', 'imp');
        catch ME
            fprintf('  RF unavailable (%s); threshold stats only.\n', ME.message);
        end
    catch ME
        fprintf('  ML classifier failed: %s\n', ME.message);
    end
end

function MultiGNSS_Dashboard(allResults, OUTPUT_ROOT)
    n = numel(allResults);
    st = cellfun(@(s) s.station, allResults, 'UniformOutput', false);
    figure('Name', 'Multi-GNSS Dashboard', 'Position', [100 50 1500 1200], 'Visible', 'off');
    % VTEC medians
    subplot(3, 2, 1); hold on; grid on;
    for s = 1:n
        try
            v = nanmedian(allResults{s}.TEC.vertical, 2);
            plot((0:numel(v)-1)/3600, v, 'LineWidth', 1.2);
        catch
        end
    end
    xlim([0 24]); xlabel('Time (UTC)'); ylabel('VTEC (TECU)');
    title('VTEC median per station'); legend(st, 'Location', 'best');
    % ROTI medians (per-second rows -> hours)
    subplot(3, 2, 2); hold on; grid on;
    for s = 1:n
        try
            r = nanmedian(allResults{s}.ROTI, 2);
            plot((0:numel(r)-1)/3600, r, 'LineWidth', 1);
            hold on;
        end
    end
    yline(0.2, 'b--', 'Minor'); yline(0.5, 'r--', 'Moderate'); yline(1.0, 'k-', 'Severe');
    xlim([0 24]); xlabel('Time (UTC)'); ylabel('ROTI (TECU/min)');
    title('ROTI median per station'); legend(st, 'Location', 'best');
    % S4 medians (amplitude scintillation) — new panel
    subplot(3, 2, 3); hold on; grid on;
    hasS4 = false;
    for s = 1:n
        try
            if isfield(allResults{s}, 'S4') && any(isfinite(allResults{s}.S4(:)))
                s4v = nanmedian(allResults{s}.S4, 2);
                plot((0:numel(s4v)-1)/3600, s4v, 'LineWidth', 1);
                hasS4 = true;
            end
        catch
        end
    end
    yline(0.15, 'b--', 'Weak');     % s4Weak
    yline(0.30, 'r--', 'Moderate'); % s4Moderate
    yline(0.50, 'k-',  'Strong');   % s4Strong
    xlim([0 24]); xlabel('Time (UTC)'); ylabel('S4 (dim-less)');
    title('S4 median per station (amplitude scintillation)');
    if hasS4, legend(st, 'Location', 'best'); end
    % S4 vs ROTI correlation (scatter, all stations combined)
    subplot(3, 2, 4); hold on; grid on;
    for s = 1:n
        try
            if ~isfield(allResults{s}, 'S4'), continue; end
            s4v  = allResults{s}.S4(:);
            rotv = allResults{s}.ROTI(:);
            ok   = isfinite(s4v) & isfinite(rotv);
            if any(ok)
                scatter(rotv(ok), s4v(ok), 2, 'filled', 'MarkerFaceAlpha', 0.3);
            end
        catch
        end
    end
    xline(0.5, 'r--', 'HandleVisibility', 'off');
    yline(0.3, 'r--', 'HandleVisibility', 'off');
    xlabel('ROTI (TECU/min)'); ylabel('S4 (dim-less)');
    title('S4 vs ROTI scatter (all stations)'); grid on;
    % Receiver DCB
    subplot(3, 2, 5);
    rcv = nan(n, 1);
    for s = 1:n
        try, rcv(s) = allResults{s}.DCB.rcv; catch, end
    end
    bar(rcv); set(gca, 'XTickLabel', st, 'XTickLabelRotation', 45);
    ylabel('TECU'); title('Receiver DCB'); grid on;
    % Station map
    subplot(3, 2, 6); hold on; grid on;
    try
        sp = loadStationPositions(fullfile(repoRoot(), 'PPPindex.txt'));
        if isempty(sp), sp = getDefaultStationData(); end
        for s = 1:n
            ix = find(strcmp({sp.name}, st{s}), 1);
            if isempty(ix), ix = find(strncmp({sp.name}, st{s}, 4), 1); end
            if ~isempty(ix)
                ll = ecef2lla(sp(ix).xyz);
                plot(ll(2), ll(1), 'k*', 'MarkerSize', 12);
                text(ll(2)+0.2, ll(1), st{s}, 'FontSize', 8);
            end
        end
    catch
    end
    xlabel('Longitude'); ylabel('Latitude'); title('Stations');
    sgtitle('Multi-GNSS TEC Dashboard (G,R,C,E,J) + S4', 'FontSize', 13);
    saveas(gcf, fullfile(OUTPUT_ROOT, 'MultiGNSS_Dashboard.png'));
    close;
end
