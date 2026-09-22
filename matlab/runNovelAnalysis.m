function Results = runNovelAnalysis(dayRoots, stations, outDir, varargin)
%{
================================================================================
 runNovelAnalysis — Master Driver for Novel Research Analysis Suite
 Q1 GPS Solutions Submission — Indian EIA Multi-GNSS Ionospheric Study
================================================================================

PURPOSE
-------
Orchestrates all 5 novel analysis modules in sequence, handles errors
gracefully, aggregates outputs into a unified Results struct, and
produces a comprehensive summary report (runNovelAnalysis_report.txt).

MODULES EXECUTED
----------------
  N1: novelAnalysis_EPBPropagation   — EPB zonal velocity field, E×B drift,
                                       bootstrap CI, amplitude-gated velocity,
                                       sunset-onset repeatability index
  N2: novelAnalysis_ScintillationCoupling — S4-ROTI coupling plane (Q1-Q4),
                                       decorrelation timescale τ_dec, MCDI,
                                       dual-frequency S4 ratio → spectral index p
  N3: novelAnalysis_TIDSpectral      — Multi-taper PSD (Thomson 1982),
                                       TID phase-speed from cross-spectrum,
                                       EIA diurnal breathing cycle analysis
  N4: novelAnalysis_ConstelDiversity — IPP Geometry Diversity Index,
                                       constellation detection efficiency η,
                                       frequency-diversity FDD and spectral index
  N5: novelAnalysis_SolarMax         — Slab thickness τ, NAI, EIA crest
                                       latitude tracking, ROTI climatology,
                                       Fejer RT-onset comparison

USAGE
-----
  % Minimal (all defaults):
  Results = runNovelAnalysis(dayRoots, stations);

  % With custom output directory and options:
  Results = runNovelAnalysis(dayRoots, stations, '/path/to/output', ...
      'runModules', [1 2 3 4 5], ...
      'IST', 5.5, ...
      'F107', [155 158 162 170 165], ...
      'Kp',   [1.5 1.8 1.2 2.3 1.0], ...
      'nBoot', 1000, ...
      'rngSeed', 42);

INPUTS
------
  dayRoots  : 1×D cell of day-root directories (same as MultiGNSS_Main)
  stations  : 1×S cell of station codes
  outDir    : [optional] output directory; default: dayRoots{1}/../ArticleFigs
  varargin  : name-value pairs (passed through to all modules):
              'runModules'  which modules to run (default [1 2 3 4 5])
              'IST'         UTC+5:30 offset (default 5.5)
              'F107'        F10.7 flux per day (default [155 158 162 170 165])
              'Kp'          Kp index per day (default [1.5 1.8 1.2 2.3 1.0])
              'nBoot'       bootstrap resamples (default 500)
              'rngSeed'     RNG seed for reproducibility (default 42)
              'elevMask'    elevation mask degrees (default 30)
              's4Mod'       S4 moderate threshold (default 0.30)
              'rotiMod'     ROTI moderate threshold (default 0.50)
              'sysList'     constellation list (default {'G','E','C','J'})
              'tidBand'     TID period band in minutes (default [15 90])
              'verbose'     print progress (default true)

OUTPUTS
-------
  Results struct:
    .N1    EPBPropagation results
    .N2    ScintillationCoupling results
    .N3    TIDSpectral results
    .N4    ConstelDiversity results
    .N5    SolarMax results
    .summary  text summary of all key findings
    .reportPath  path to written report file
    .timing   per-module wall-clock times (seconds)
    .errors   cell array of error messages (empty if all OK)

REQUIRED FILES (all in same directory as this script)
------
  novelAnalysis_EPBPropagation.m
  novelAnalysis_ScintillationCoupling.m
  novelAnalysis_TIDSpectral.m
  novelAnalysis_ConstelDiversity.m
  novelAnalysis_SolarMax.m

EXAMPLE (integration into MultiGNSS_Main.m)
--------------------------------------------
  % Add at end of MultiGNSS_Main pipeline:
  if cfg.runNovelAnalysis
      Results = runNovelAnalysis(dayRoots, cfg.stations, cfg.outRoot, ...
          'F107', cfg.F107, 'Kp', cfg.Kp);
  end

Q1 COMPLIANCE
-------------
  - All 5 modules produce article-ready 300 DPI 4-panel figures
  - All outputs include .mat + .csv for reproducibility
  - RNG seed 42 fixed across all bootstrap operations
  - No MATLAB Statistics/Signal Processing Toolbox required
  - Full error trapping: one failing module does not abort others

================================================================================
%}

    tStart = tic;

    %% -------- defaults and parse ------------------------------------------
    p.runModules = [1 2 3 4 5];
    p.IST        = 5.5;
    p.F107       = [155, 158, 162, 170, 165];
    p.Kp         = [1.5,  1.8,  1.2,  2.3,  1.0];
    p.nBoot      = 500;
    p.rngSeed    = 42;
    p.elevMask   = 30;
    p.s4Mod      = 0.30;
    p.rotiMod    = 0.50;
    p.sysList    = {'G','E','C','J'};
    p.tidBand    = [15 90];
    p.verbose    = true;

    for k = 1:2:numel(varargin)
        if isfield(p, varargin{k}), p.(varargin{k}) = varargin{k+1}; end
    end

    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % Common options to pass to all modules
    commonOpts = {'IST', p.IST, 'nBoot', p.nBoot, 'rngSeed', p.rngSeed, ...
        'elevMask', p.elevMask, 's4Mod', p.s4Mod, 'rotiMod', p.rotiMod};

    if p.verbose
        fprintf('\n');
        fprintf('╔══════════════════════════════════════════════════════════╗\n');
        fprintf('║  runNovelAnalysis — Q1 Novel Research Analysis Suite     ║\n');
        fprintf('║  Indian EIA Multi-GNSS Ionospheric Study                 ║\n');
        fprintf('║  Modules: %s                              ║\n', ...
            num2str(p.runModules, '%d '));
        fprintf('╚══════════════════════════════════════════════════════════╝\n');
        fprintf('  Day roots: %d days\n', numel(dayRoots));
        fprintf('  Stations : %s\n', strjoin(stations, ', '));
        fprintf('  Output   : %s\n\n', outDir);
    end

    Results = struct('N1',[],'N2',[],'N3',[],'N4',[],'N5',[], ...
        'summary','','reportPath','','timing',[],'errors',{{}});
    timing = nan(1,5);
    errors = {};

    %% ================================================================
    %% MODULE N1: EPB Propagation
    %% ================================================================
    if ismember(1, p.runModules)
        if p.verbose, fprintf('━━━ N1: EPB Propagation (velocity, E×B drift, repeatability) ━━━\n'); end
        t1 = tic;
        try
            Results.N1 = novelAnalysis_EPBPropagation(dayRoots, stations, outDir, ...
                commonOpts{:});
            timing(1) = toc(t1);
            if p.verbose
                fprintf('[N1] ✓ Completed in %.1f s | Figure: %s\n\n', timing(1), Results.N1.figPath);
            end
        catch ME
            timing(1) = toc(t1);
            msg = sprintf('[N1] ERROR: %s (%s line %d)', ME.message, ME.stack(1).name, ME.stack(1).line);
            fprintf('%s\n', msg);
            errors{end+1} = msg;
        end
    end

    %% ================================================================
    %% MODULE N2: Scintillation Coupling
    %% ================================================================
    if ismember(2, p.runModules)
        if p.verbose, fprintf('━━━ N2: Scintillation Coupling (S4-ROTI quadrant, MCDI, τ_dec) ━━━\n'); end
        t1 = tic;
        try
            Results.N2 = novelAnalysis_ScintillationCoupling(dayRoots, stations, outDir, ...
                commonOpts{:});
            timing(2) = toc(t1);
            if p.verbose
                fprintf('[N2] ✓ Completed in %.1f s | Figure: %s\n\n', timing(2), Results.N2.figPath);
            end
        catch ME
            timing(2) = toc(t1);
            msg = sprintf('[N2] ERROR: %s (%s line %d)', ME.message, ME.stack(1).name, ME.stack(1).line);
            fprintf('%s\n', msg);
            errors{end+1} = msg;
        end
    end

    %% ================================================================
    %% MODULE N3: TID Spectral Analysis
    %% ================================================================
    if ismember(3, p.runModules)
        if p.verbose, fprintf('━━━ N3: TID Spectral (multi-taper PSD, phase-speed, EIA breathing) ━━━\n'); end
        t1 = tic;
        try
            Results.N3 = novelAnalysis_TIDSpectral(dayRoots, stations, outDir, ...
                commonOpts{:}, 'tidBand', p.tidBand);
            timing(3) = toc(t1);
            if p.verbose
                fprintf('[N3] ✓ Completed in %.1f s | Figure: %s\n\n', timing(3), Results.N3.figPath);
            end
        catch ME
            timing(3) = toc(t1);
            msg = sprintf('[N3] ERROR: %s (%s line %d)', ME.message, ME.stack(1).name, ME.stack(1).line);
            fprintf('%s\n', msg);
            errors{end+1} = msg;
        end
    end

    %% ================================================================
    %% MODULE N4: Constellation Diversity
    %% ================================================================
    if ismember(4, p.runModules)
        if p.verbose, fprintf('━━━ N4: Constellation Diversity (GDI, detection η, freq-diversity) ━━━\n'); end
        t1 = tic;
        try
            Results.N4 = novelAnalysis_ConstelDiversity(dayRoots, stations, outDir, ...
                commonOpts{:}, 'sysList', p.sysList);
            timing(4) = toc(t1);
            if p.verbose
                fprintf('[N4] ✓ Completed in %.1f s | Figure: %s\n\n', timing(4), Results.N4.figPath);
            end
        catch ME
            timing(4) = toc(t1);
            msg = sprintf('[N4] ERROR: %s (%s line %d)', ME.message, ME.stack(1).name, ME.stack(1).line);
            fprintf('%s\n', msg);
            errors{end+1} = msg;
        end
    end

    %% ================================================================
    %% MODULE N5: Solar Maximum Characterization
    %% ================================================================
    if ismember(5, p.runModules)
        if p.verbose, fprintf('━━━ N5: Solar Maximum (slab thickness, NAI, EIA crest, ROTI clim) ━━━\n'); end
        t1 = tic;
        try
            Results.N5 = novelAnalysis_SolarMax(dayRoots, stations, outDir, ...
                commonOpts{:}, 'F107', p.F107, 'Kp', p.Kp);
            timing(5) = toc(t1);
            if p.verbose
                fprintf('[N5] ✓ Completed in %.1f s | Figure: %s\n\n', timing(5), Results.N5.figPath);
            end
        catch ME
            timing(5) = toc(t1);
            msg = sprintf('[N5] ERROR: %s (%s line %d)', ME.message, ME.stack(1).name, ME.stack(1).line);
            fprintf('%s\n', msg);
            errors{end+1} = msg;
        end
    end

    totalTime = toc(tStart);
    Results.timing = timing;
    Results.errors = errors;

    %% ================================================================
    %% COMPILE SUMMARY REPORT
    %% ================================================================
    sumLines = {};
    sumLines{end+1} = '================================================================================';
    sumLines{end+1} = ' NOVEL ANALYSIS SUITE — RESULTS SUMMARY';
    sumLines{end+1} = ' Indian EIA Multi-GNSS Study, 1-5 May 2025 (Solar Maximum)';
    sumLines{end+1} = '================================================================================';
    sumLines{end+1} = sprintf(' Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    sumLines{end+1} = sprintf(' Total runtime: %.1f s', totalTime);
    sumLines{end+1} = '';

    % N1 summary
    if ~isempty(Results.N1)
        sumLines{end+1} = '--- N1: EPB Zonal Propagation ---';
        try
            if isfield(Results.N1, 'velMed')
                sumLines{end+1} = sprintf('  Median EPB zonal velocity  : %.0f m/s', Results.N1.velMed*1000);
            end
            if isfield(Results.N1, 'velCI95')
                sumLines{end+1} = sprintf('  Bootstrap 95%% CI           : [%.0f, %.0f] m/s', ...
                    Results.N1.velCI95(1)*1000, Results.N1.velCI95(2)*1000);
            end
            if isfield(Results.N1, 'repeatabilitySTD')
                sumLines{end+1} = sprintf('  Sunset-onset STD           : %.2f h', Results.N1.repeatabilitySTD);
            end
        catch, end
        sumLines{end+1} = '';
    end

    % N2 summary
    if ~isempty(Results.N2)
        sumLines{end+1} = '--- N2: S4-ROTI Scintillation Coupling ---';
        try
            if isfield(Results.N2, 'quadrantFrac')
                qf = Results.N2.quadrantFrac;
                sumLines{end+1} = sprintf('  Quadrant fractions Q1/Q2/Q3/Q4: %.1f%%/%.1f%%/%.1f%%/%.1f%%', ...
                    qf(1)*100, qf(2)*100, qf(3)*100, qf(4)*100);
            end
            if isfield(Results.N2, 'tau_dec_EPB') && isfield(Results.N2,'tau_dec_quiet')
                sumLines{end+1} = sprintf('  τ_dec: EPB=%.1f s  Quiet=%.1f s', ...
                    Results.N2.tau_dec_EPB, Results.N2.tau_dec_quiet);
            end
            if isfield(Results.N2, 'MCDI_max')
                sumLines{end+1} = sprintf('  MCDI peak (4 May EPB)      : %.3f', Results.N2.MCDI_max);
            end
        catch, end
        sumLines{end+1} = '';
    end

    % N3 summary
    if ~isempty(Results.N3)
        sumLines{end+1} = '--- N3: TID Spectral Analysis ---';
        try
            if isfield(Results.N3, 'tidCatalog') && numel(Results.N3.tidCatalog) > 0
                sumLines{end+1} = sprintf('  TID events detected        : %d', numel(Results.N3.tidCatalog));
            end
            if isfield(Results.N3,'phaseSpeed') && isfield(Results.N3.phaseSpeed,'speedMed')
                sumLines{end+1} = sprintf('  Median TID phase speed     : %.0f m/s', ...
                    Results.N3.phaseSpeed.speedMed*1000);
            end
            if isfield(Results.N3,'eiaBreath')
                eb = Results.N3.eiaBreath;
                Amps = [eb.A_TECU]; ok = logical([eb.fitOK]);
                if any(ok)
                    sumLines{end+1} = sprintf('  EIA breathing amplitude    : %.1f ± %.1f TECU', ...
                        nanmean(Amps(ok)), nanstd(Amps(ok)));
                end
            end
        catch, end
        sumLines{end+1} = '';
    end

    % N4 summary
    if ~isempty(Results.N4)
        sumLines{end+1} = '--- N4: Constellation Diversity ---';
        try
            if isfield(Results.N4,'detection') && ~isempty(Results.N4.detection)
                de = Results.N4.detection;
                for k = 1:numel(de)
                    sumLines{end+1} = sprintf('  %s detection η=%.3f  FA=%.3f', ...
                        de(k).sys, de(k).eta, de(k).falseAlarm);
                end
            end
            if isfield(Results.N4,'freqDiv') && isfield(Results.N4.freqDiv,'welchT_GvsE')
                sumLines{end+1} = sprintf('  Spectral index G vs E: Welch t=%.2f p=%.4f', ...
                    Results.N4.freqDiv.welchT_GvsE, Results.N4.freqDiv.welchP_GvsE);
            end
        catch, end
        sumLines{end+1} = '';
    end

    % N5 summary
    if ~isempty(Results.N5)
        sumLines{end+1} = '--- N5: Solar Maximum Characterization ---';
        try
            if isfield(Results.N5,'rotiClim') && isfield(Results.N5.rotiClim,'fejerOnset')
                fo = Results.N5.rotiClim.fejerOnset;
                oo = Results.N5.rotiClim.obsEPBonset;
                drt= Results.N5.rotiClim.deltaRT;
                sumLines{end+1} = sprintf('  Mean Fejer RT onset        : %.1f h LT', nanmean(fo));
                sumLines{end+1} = sprintf('  Mean observed EPB onset    : %.1f h LT', ...
                    nanmean(oo(isfinite(oo))));
                sumLines{end+1} = sprintf('  Mean ΔRT (obs-pred)        : %.1f h', ...
                    nanmean(drt(isfinite(drt))));
            end
            if isfield(Results.N5,'crestDrift')
                sumLines{end+1} = sprintf('  EIA crest drift rate       : %.3f ± %.3f deg/h', ...
                    nanmean(Results.N5.crestDrift), nanstd(Results.N5.crestDrift));
            end
        catch, end
        sumLines{end+1} = '';
    end

    % Errors
    if ~isempty(errors)
        sumLines{end+1} = '--- ERRORS ---';
        for k = 1:numel(errors)
            sumLines{end+1} = ['  ' errors{k}];
        end
        sumLines{end+1} = '';
    end

    sumLines{end+1} = '--- MODULE TIMING ---';
    moduleNames = {'N1:EPBPropagation','N2:ScintCoupling','N3:TIDSpectral','N4:ConstelDiversity','N5:SolarMax'};
    for k = 1:5
        if isfinite(timing(k))
            sumLines{end+1} = sprintf('  %-22s : %.1f s', moduleNames{k}, timing(k));
        end
    end
    sumLines{end+1} = sprintf('  %-22s : %.1f s', 'TOTAL', totalTime);
    sumLines{end+1} = '';
    sumLines{end+1} = '--- OUTPUT FILES ---';
    for mod = {Results.N1, Results.N2, Results.N3, Results.N4, Results.N5}
        m = mod{1};
        if ~isempty(m)
            if isfield(m,'figPath') && ~isempty(m.figPath)
                sumLines{end+1} = sprintf('  Figure : %s', m.figPath);
            end
            if isfield(m,'matPath') && ~isempty(m.matPath)
                sumLines{end+1} = sprintf('  MAT    : %s', m.matPath);
            end
            if isfield(m,'csvPath') && ~isempty(m.csvPath)
                sumLines{end+1} = sprintf('  CSV    : %s', m.csvPath);
            end
        end
    end
    sumLines{end+1} = '================================================================================';

    sumText = strjoin(sumLines, '\n');
    Results.summary = sumText;

    % Write report
    reportPath = fullfile(outDir, 'runNovelAnalysis_report.txt');
    try
        fid = fopen(reportPath, 'w');
        fprintf(fid, '%s\n', sumText);
        fclose(fid);
        Results.reportPath = reportPath;
    catch
        Results.reportPath = '';
    end

    %% ================================================================
    %% FINAL PRINT
    %% ================================================================
    if p.verbose
        fprintf('\n');
        fprintf('╔══════════════════════════════════════════════════════════╗\n');
        fprintf('║  runNovelAnalysis — COMPLETE                             ║\n');
        fprintf('╠══════════════════════════════════════════════════════════╣\n');
        fprintf('║  Total runtime : %-6.1f s                                ║\n', totalTime);
        fprintf('║  Modules OK    : %-2d / %-2d                                ║\n', ...
            sum(isfinite(timing)), sum(ismember([1:5], p.runModules)));
        fprintf('║  Errors        : %-2d                                     ║\n', numel(errors));
        fprintf('║  Report        : %s\n', reportPath);
        fprintf('╚══════════════════════════════════════════════════════════╝\n\n');

        if ~isempty(errors)
            fprintf('[!] The following modules encountered errors:\n');
            for k = 1:numel(errors)
                fprintf('    %s\n', errors{k});
            end
            fprintf('\n');
        end
    end
end
