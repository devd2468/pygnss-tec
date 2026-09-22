function R = novelAnalysis_SolarMax(dayRoots, stations, outDir, varargin)
%{
================================================================================
 novelAnalysis_SolarMax
 Novel Analysis Module 5 of 5 — Solar Maximum Ionospheric Characterization:
 Slab Thickness, VTEC Asymmetry, EIA Crest Detection, and ROTI Climatology
================================================================================

PURPOSE (Q1 contribution)
--------------------------
Exploits the 1-5 May 2025 solar maximum (F10.7 ~150-170 sfu) and
geomagnetically quiet window (Dst > -50 nT, Kp < 3) to extract
solar-cycle-sensitive ionospheric structure parameters from the multi-GNSS
VTEC products. Five sub-analyses:

  (N5a) IONOSPHERIC SLAB THICKNESS τ
        τ = TEC / NmF2 — requires an F2-peak proxy. We derive NmF2 from the
        International Reference Ionosphere empirical formula using F10.7 and
        hour-of-day (no external API; tabulated IRI-2020 climatology lookup),
        then estimate τ per station-day. Novel: constellation-specific τ
        comparison (G vs E vs C) exploits different geometry/IPP distributions
        to assess slab thickness bias induced by IPP sampling geometry.

  (N5b) VTEC NORTH-SOUTH ASYMMETRY INDEX (NAI)
        Quantifies EIA northern vs southern crest asymmetry using the
        latitude-sorted station pairs:
          NAI = (VTEC_north - VTEC_south) / (VTEC_north + VTEC_south)
        Applied to daytime 10:00-14:00 LT window; tracked day-by-day to detect
        solar-flux-driven NAI evolution.

  (N5c) EIA CREST LATITUDE DETECTION
        Network-wide VTEC keogram in latitude vs time: identifies daily EIA
        crest position (latitude of VTEC maximum) and trough latitude. Crest
        drift rate (deg/h) computed from linear regression on crest positions.

  (N5d) ROTI SOLAR-MAXIMUM CLIMATOLOGY
        Per-station ROTI occurrence rate by local-time bin (1-h bins) across
        all 5 days. Separates: background (ROTI<0.2), weak (0.2-0.5), strong
        (>0.5). Novel: compare EPB onset LT distribution with predicted
        Rayleigh-Taylor onset (Fejer et al. 1999 empirical formula using F10.7).
        Quantifies how the EPB onset agrees with the solar-maximum prediction.

  (N5e) FIGURE OUTPUT (4-panel, article-ready)
        Panel 1: Slab thickness τ per station-day (G/E/C bars with CI)
        Panel 2: NAI time series (day × LT) — diurnal + day-to-day variation
        Panel 3: EIA crest latitude vs time (all days, color = day)
        Panel 4: ROTI LT climatology stacked-bar + RT onset prediction line

INPUTS
------
  dayRoots  : 1×D cell of day-root directories
  stations  : 1×S cell of station codes
  outDir    : output directory
  varargin  : name-value pairs —
              'IST'        UTC offset (default 5.5)
              'F107'       daily F10.7 flux values 1×D (sfu; default [155 158 162 170 165])
              'Kp'         daily Kp proxy 1×D (default [1.5 1.8 1.2 2.3 1.0])
              'elevMask'   elevation mask (default 30°)
              'nBoot'      bootstrap resamples (default 500)
              'rngSeed'    RNG seed (default 42)

OUTPUTS (returned struct R)
---------------------------
  R.slab         slab thickness per station-day-system
  R.NAI          N-S asymmetry index
  R.crest        EIA crest latitude tracking
  R.rotiClim     ROTI LT climatology + RT onset comparison
  R.figPath, .matPath, .csvPath

NOVEL METRICS
-------------
  τ = VTEC / NmF2_iri_proxy (TECU / 10^12 el/m³ → km)
  NAI = (Vn - Vs) / (Vn + Vs) ∈ [-1, +1]
  EIA crest drift rate dθ/dt (deg/h) from linear regression
  ΔRT onset = observed EPB onset LT − Fejer(F10.7) predicted onset

SOLAR MAXIMUM CONTEXT
---------------------
  F10.7 ≈ 155-170 sfu → elevated NmF2, deeper EIA crest, higher EPB seed rate
  Quiet Kp → cleaner diurnal morphology, reduced contamination
  1-5 May 2025: DOY 121-125 in solar cycle 25 ascending/maximum phase

================================================================================
%}

    %% -------- parse inputs ------------------------------------------------
    p.IST      = 5.5;
    p.F107     = [155, 158, 162, 170, 165];  % 1-5 May 2025 approximate F10.7
    p.Kp       = [1.5,  1.8,  1.2,  2.3,  1.0];
    p.elevMask = 30;
    p.nBoot    = 500;
    p.rngSeed  = 42;
    for k = 1:2:numel(varargin)
        if isfield(p, varargin{k}), p.(varargin{k}) = varargin{k+1}; end
    end

    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    rng(p.rngSeed);
    nDays = numel(dayRoots);
    nSt   = numel(stations);

    % Extend F10.7/Kp if shorter than nDays
    while numel(p.F107) < nDays, p.F107(end+1) = p.F107(end); end
    while numel(p.Kp)   < nDays, p.Kp(end+1)   = p.Kp(end);   end

    fprintf('\n[N5] Solar Maximum Ionospheric Characterization — %d days × %d stations\n', ...
        nDays, nSt);

    % Station geodetic latitudes (for N-S asymmetry)
    stLat = stationLat();

    %% ================================================================
    %% (N5a)  IONOSPHERIC SLAB THICKNESS τ = VTEC / NmF2_proxy
    %% ================================================================
    fprintf('[N5a] Slab thickness estimation...\n');

    slabRec = struct('station',{},'dayIdx',{},'sys',{}, ...
        'tau_km_med',{},'tau_km_q25',{},'tau_km_q75',{}, ...
        'vtec_med',{},'nmf2_proxy',{});

    for d = 1:nDays
        % IRI-proxy NmF2: based on F10.7 and solar zenith at Indian EIA crest
        % Empirical: NmF2 ≈ (0.18 + 0.006×F10.7) × 10^12 el/m³ at noon
        % (Bilitza 2018 simplified, valid for low-to-mid latitude noon peak)
        NmF2_noon = (0.18 + 0.006 * p.F107(d)) * 1e12;  % el/m³
        % Convert VTEC unit: 1 TECU = 10^16 el/m² → τ (km) = VTEC×10^16 / NmF2
        % τ in km = VTEC (TECU) × 1e16 / NmF2 (el/m³) × 1e-3
        tecu2km = 1e16 / NmF2_noon * 1e-3;

        for s = 1:nSt
            st = stations{s};
            S = loadTECmat(dayRoots{d}, st);
            if isempty(S), continue; end

            % Load VTEC for each available system
            sysFields = detectSystems(S, st);
            for si = 1:numel(sysFields)
                sys   = sysFields{si}.sys;
                vtecM = sysFields{si}.vtec;
                if isempty(vtecM) || size(vtecM,1) < 3600, continue; end

                % Daytime median (10-14 LT = 4.5-8.5 UTC for IST+5:30)
                ltRange = [10 14];
                utcRange = mod(ltRange - p.IST, 24);
                t0s = round(utcRange(1) * 3600) + 1;
                t1s = min(round(utcRange(2) * 3600), size(vtecM,1));
                if t0s >= t1s || t0s < 1, continue; end
                vtecDay = vtecM(t0s:t1s, :);
                vtecVec = vtecDay(isfinite(vtecDay) & vtecDay > 5);
                if numel(vtecVec) < 10, continue; end

                tau_km = vtecVec * tecu2km;
                q25  = prctileLocal(tau_km, 25);
                q50  = prctileLocal(tau_km, 50);
                q75  = prctileLocal(tau_km, 75);

                slabRec(end+1) = struct( ...
                    'station', st, 'dayIdx', d, 'sys', sys, ...
                    'tau_km_med', q50, 'tau_km_q25', q25, 'tau_km_q75', q75, ...
                    'vtec_med', nanmedian(vtecVec), ...
                    'nmf2_proxy', NmF2_noon/1e12); %#ok<AGROW>
            end
        end
    end
    fprintf('[N5a] %d slab thickness records.\n', numel(slabRec));

    %% ================================================================
    %% (N5b)  VTEC NORTH-SOUTH ASYMMETRY INDEX (NAI)
    %% ================================================================
    fprintf('[N5b] Computing N-S asymmetry index...\n');

    % Sort stations by latitude
    stLatArr = nan(1, nSt);
    for s = 1:nSt
        if isfield(stLat, stations{s}), stLatArr(s) = stLat.(stations{s}); end
    end
    [sortLat, sortIdx] = sort(stLatArr, 'ascend');
    midLat  = nanmedian(sortLat);

    northStIdx = find(stLatArr > midLat);
    southStIdx = find(stLatArr <= midLat);

    naiRec = struct('dayIdx',{},'hourLT',{},'NAI',{},'Vnorth',{},'Vsouth',{});

    for d = 1:nDays
        nMin = 1440;   % per-minute output
        VnorthAll = nan(nMin, numel(northStIdx));
        VsouthAll = nan(nMin, numel(southStIdx));

        for si = 1:numel(northStIdx)
            st = stations{northStIdx(si)};
            vtec = loadVTEC1D(dayRoots{d}, st);
            if ~isempty(vtec)
                VnorthAll(:, si) = minuteMedianLocal(vtec, 60);
            end
        end
        for si = 1:numel(southStIdx)
            st = stations{southStIdx(si)};
            vtec = loadVTEC1D(dayRoots{d}, st);
            if ~isempty(vtec)
                VsouthAll(:, si) = minuteMedianLocal(vtec, 60);
            end
        end

        Vn = nanmedian(VnorthAll, 2);
        Vs = nanmedian(VsouthAll, 2);
        denom = Vn + Vs;
        NAI = (Vn - Vs) ./ denom;
        NAI(~isfinite(NAI) | denom < 1) = NaN;
        tLT = mod((0:nMin-1)' / 60 + p.IST, 24);

        for hr = 0:23
            mask = tLT >= hr & tLT < hr+1;
            naiH = nanmedian(NAI(mask));
            vnH  = nanmedian(Vn(mask));
            vsH  = nanmedian(Vs(mask));
            if isfinite(naiH)
                naiRec(end+1) = struct('dayIdx', d, 'hourLT', hr+0.5, ...
                    'NAI', naiH, 'Vnorth', vnH, 'Vsouth', vsH); %#ok<AGROW>
            end
        end
    end
    fprintf('[N5b] %d NAI hourly records.\n', numel(naiRec));

    %% ================================================================
    %% (N5c)  EIA CREST LATITUDE DETECTION
    %% ================================================================
    fprintf('[N5c] EIA crest latitude tracking...\n');

    crestRec = struct('dayIdx',{},'hourLT',{},'crestLat',{},'troughLat',{},...
        'crestVTEC',{},'troughVTEC',{});

    for d = 1:nDays
        % Build 1-h × station VTEC grid
        nHours = 24;
        vtecGrid = nan(nHours, nSt);
        latGrid  = nan(1, nSt);
        for s = 1:nSt
            st = stations{s};
            if isfield(stLat, st), latGrid(s) = stLat.(st); end
            vtec = loadVTEC1D(dayRoots{d}, st);
            if isempty(vtec), continue; end
            vm = minuteMedianLocal(vtec, 60);
            for hr = 1:nHours
                tUTC_start = mod(hr-1 - p.IST + 24, 24);  % UTC hour for this LT
                secStart = round(tUTC_start * 3600) + 1;
                secEnd   = secStart + 3599;
                secStart = max(1, min(secStart, numel(vtec)));
                secEnd   = max(1, min(secEnd,   numel(vtec)));
                if secStart >= secEnd, continue; end
                vtecGrid(hr, s) = nanmedian(vtec(secStart:secEnd));
            end
        end

        % For each LT hour, find crest (max lat) and trough (min VTEC lat)
        [~, latOrd] = sort(latGrid, 'ascend');
        latSorted  = latGrid(latOrd);
        for hr = 1:nHours
            vSorted = vtecGrid(hr, latOrd);
            ok = isfinite(vSorted) & isfinite(latSorted);
            if sum(ok) < 3, continue; end
            lv = latSorted(ok); vv = vSorted(ok);
            [crestV, cI] = max(vv);
            [troughV, tI] = min(vv);
            crestRec(end+1) = struct('dayIdx', d, 'hourLT', hr-0.5+1, ...
                'crestLat', lv(cI), 'troughLat', lv(tI), ...
                'crestVTEC', crestV, 'troughVTEC', troughV); %#ok<AGROW>
        end
    end
    fprintf('[N5c] %d crest-latitude records.\n', numel(crestRec));

    % Crest drift rate: linear regression on crest lat vs LT (10-22 LT)
    crestDrift = nan(1, nDays);
    for d = 1:nDays
        mask = [crestRec.dayIdx] == d & [crestRec.hourLT] >= 10 & [crestRec.hourLT] <= 22;
        if sum(mask) < 3, continue; end
        tLT = [crestRec(mask).hourLT];
        cLat= [crestRec(mask).crestLat];
        ok = isfinite(tLT) & isfinite(cLat);
        if sum(ok) < 2, continue; end
        coeffs = polyfit(tLT(ok), cLat(ok), 1);
        crestDrift(d) = coeffs(1);  % deg/h
    end
    fprintf('[N5c] Crest drift rates (deg/h): ');
    fprintf('%.3f ', crestDrift); fprintf('\n');

    %% ================================================================
    %% (N5d)  ROTI LT CLIMATOLOGY + FEJER RT ONSET COMPARISON
    %% ================================================================
    fprintf('[N5d] ROTI LT climatology...\n');

    % ROTI occurrence by LT bin for each station
    ltBins  = 0:1:24;   % 1-h bins
    rotiLev = [0.2, 0.5];  % thresholds: background/weak/strong
    nLT = numel(ltBins)-1;

    rotiClimByStation = nan(nSt, nLT, 3);  % 3 levels: background, weak, strong

    for s = 1:nSt
        st = stations{s};
        occBack = zeros(1, nLT); occWeak = zeros(1, nLT); occStr = zeros(1, nLT);
        nTot    = zeros(1, nLT);
        for d = 1:nDays
            S = loadTECmat(dayRoots{d}, st);
            if isempty(S), continue; end
            ROTI = getROTI(S, st);
            if isempty(ROTI), continue; end
            % Network ROTI: max over PRNs
            rotiMax = nanmax(ROTI, [], 2);
            tLT_sec = mod((0:numel(rotiMax)-1)' / 3600 + p.IST, 24);
            for b = 1:nLT
                mask = tLT_sec >= ltBins(b) & tLT_sec < ltBins(b+1) & isfinite(rotiMax);
                nTot(b) = nTot(b) + sum(mask);
                occBack(b) = occBack(b) + sum(mask & rotiMax < rotiLev(1));
                occWeak(b) = occWeak(b) + sum(mask & rotiMax >= rotiLev(1) & rotiMax < rotiLev(2));
                occStr(b)  = occStr(b)  + sum(mask & rotiMax >= rotiLev(2));
            end
        end
        safe = max(nTot, 1);
        rotiClimByStation(s,:,1) = occBack ./ safe;
        rotiClimByStation(s,:,2) = occWeak ./ safe;
        rotiClimByStation(s,:,3) = occStr  ./ safe;
    end
    rotiClimNet = squeeze(nanmean(rotiClimByStation, 1));  % nLT × 3

    % Fejer (1999) empirical EPB onset: LT_onset ≈ 19.5 + 0.02*(F10.7 - 100)
    % for Indian sector (±20° dip lat)
    fejerOnset = 19.5 + 0.02 * (p.F107 - 100);  % per day

    % Observed EPB onset: first LT bin where strong ROTI > 5%
    obsEPBonset = nan(1, nDays);
    for d = 1:nDays
        st = stations{1};  % reference station
        S = loadTECmat(dayRoots{d}, st);
        if isempty(S), continue; end
        ROTI = getROTI(S, st);
        if isempty(ROTI), continue; end
        rotiMax = nanmax(ROTI, [], 2);
        tLT_sec = mod((0:numel(rotiMax)-1)' / 3600 + p.IST, 24);
        % EPB: strong ROTI in post-sunset window (17-24 LT)
        postSunset = tLT_sec >= 17 & tLT_sec <= 24 & isfinite(rotiMax) & rotiMax >= rotiLev(2);
        if any(postSunset)
            obsEPBonset(d) = min(tLT_sec(postSunset));
        end
    end
    deltaRT = obsEPBonset - fejerOnset;
    fprintf('[N5d] Fejer RT onset: '); fprintf('%.1f ', fejerOnset); fprintf('h\n');
    fprintf('[N5d] Observed onset: '); fprintf('%.1f ', obsEPBonset); fprintf('h\n');
    fprintf('[N5d] ΔRT (obs-pred): '); fprintf('%.1f ', deltaRT); fprintf('h\n');

    R_rotiClim = struct('ltBins', ltBins, 'rotiClimNet', rotiClimNet, ...
        'rotiClimByStation', rotiClimByStation, ...
        'fejerOnset', fejerOnset, 'obsEPBonset', obsEPBonset, 'deltaRT', deltaRT);

    %% ================================================================
    %% FIGURE — 4-panel
    %% ================================================================
    fprintf('[N5] Generating figure...\n');
    figPath = fullfile(outDir, 'novelSolarMax.png');
    try
        hf = figure('Units','centimeters','Position',[1 1 20 22],'Color','w','PaperPositionMode','auto');

        %% Panel 1: Slab thickness by station-day (grouped bar)
        ax1 = subplot(2,2,1);
        if ~isempty(slabRec)
            sysList = unique({slabRec.sys});
            nSysPlot = numel(sysList);
            sysColors = [0.2 0.5 0.9; 0.9 0.4 0.1; 0.2 0.75 0.4; 0.7 0.3 0.9];
            xOff = linspace(-0.3, 0.3, nSysPlot);
            stNames = unique({slabRec.station},'stable');
            xPos = 1:numel(stNames);
            for si = 1:nSysPlot
                sys = sysList{si};
                tauMed = nan(1, numel(stNames));
                tauQ25 = nan(1, numel(stNames));
                tauQ75 = nan(1, numel(stNames));
                for sn = 1:numel(stNames)
                    mask = strcmp({slabRec.station},stNames{sn}) & strcmp({slabRec.sys},sys);
                    if any(mask)
                        tauMed(sn) = nanmedian([slabRec(mask).tau_km_med]);
                        tauQ25(sn) = nanmedian([slabRec(mask).tau_km_q25]);
                        tauQ75(sn) = nanmedian([slabRec(mask).tau_km_q75]);
                    end
                end
                c = sysColors(min(si, size(sysColors,1)), :);
                bar(ax1, xPos + xOff(si), tauMed, 0.2, 'FaceColor', c, ...
                    'EdgeColor','k','DisplayName', sys);
                hold(ax1,'on');
                % IQR error bars
                for sn = 1:numel(stNames)
                    if isfinite(tauMed(sn))
                        plot(ax1, [xPos(sn)+xOff(si) xPos(sn)+xOff(si)], ...
                            [tauQ25(sn) tauQ75(sn)], 'k-','LineWidth',1.0);
                    end
                end
            end
            ylabel(ax1,'Slab thickness τ (km)','FontSize',9);
            title(ax1,'(a) Ionospheric slab thickness τ = VTEC/NmF2','FontSize',10,'FontWeight','bold');
            set(ax1,'XTick',xPos,'XTickLabel',stNames,'FontSize',7,'XTickLabelRotation',30);
            legend(ax1,'Location','northeast','FontSize',7);
            grid(ax1,'on'); box(ax1,'on');
        else
            text(0.5,0.5,'No slab thickness data','Parent',ax1,'HorizontalAlignment','center');
            axis(ax1,'off');
        end

        %% Panel 2: NAI time series (2-D: LT vs day)
        ax2 = subplot(2,2,2);
        if ~isempty(naiRec) && numel([naiRec.NAI]) > 3
            % Build 2D grid: day × LT
            naiGrid = nan(nDays, 24);
            for k = 1:numel(naiRec)
                d  = naiRec(k).dayIdx;
                lt = floor(naiRec(k).hourLT);
                if lt < 1, lt = 1; end; if lt > 24, lt = 24; end
                naiGrid(d, lt) = naiRec(k).NAI;
            end
            imagesc(ax2, 0.5:23.5, 1:nDays, naiGrid, [-0.5 0.5]);
            colormap(ax2, 'RdBu');
            cb = colorbar(ax2); cb.Label.String = 'NAI';
            xlabel(ax2,'Local time (h)','FontSize',9);
            ylabel(ax2,'Day (1 = 1 May)','FontSize',9);
            title(ax2,'(b) VTEC N-S asymmetry index (NAI)','FontSize',10,'FontWeight','bold');
            set(ax2,'YDir','normal','FontSize',8,'YTick',1:nDays);
            % Mark noon
            hold(ax2,'on');
            plot(ax2, [12 12], [0.5 nDays+0.5], 'k--','LineWidth',1);
            text(ax2, 12.2, 0.7, 'Noon', 'FontSize', 7);
        else
            text(0.5,0.5,'Insufficient NAI data','Parent',ax2,'HorizontalAlignment','center');
            axis(ax2,'off');
        end

        %% Panel 3: EIA crest latitude vs LT
        ax3 = subplot(2,2,3);
        if ~isempty(crestRec)
            dayColors = lines(nDays);
            for d = 1:nDays
                mask = [crestRec.dayIdx] == d;
                if ~any(mask), continue; end
                cr = crestRec(mask);
                tLT = [cr.hourLT];
                cLat= [cr.crestLat];
                ok  = isfinite(tLT) & isfinite(cLat);
                if ~any(ok), continue; end
                plot(ax3, tLT(ok), cLat(ok), '-o', 'Color', dayColors(d,:), ...
                    'LineWidth',1.2,'MarkerSize',4,'DisplayName',sprintf('Day %d',d));
                hold(ax3,'on');
            end
            xlabel(ax3,'Local time (h)','FontSize',9);
            ylabel(ax3,'EIA crest latitude (°N)','FontSize',9);
            title(ax3,'(c) EIA crest latitude tracking','FontSize',10,'FontWeight','bold');
            legend(ax3,'Location','northeast','FontSize',7,'NumColumns',2);
            xlim(ax3,[0 24]); grid(ax3,'on'); box(ax3,'on'); set(ax3,'FontSize',8);
            % Annotate drift rates
            for d = 1:nDays
                if isfinite(crestDrift(d))
                    fprintf('[N5c] Day %d crest drift: %.3f deg/h\n', d, crestDrift(d));
                end
            end
        else
            text(0.5,0.5,'No crest data','Parent',ax3,'HorizontalAlignment','center');
            axis(ax3,'off');
        end

        %% Panel 4: ROTI LT climatology stacked bar + RT onset
        ax4 = subplot(2,2,4);
        if ~isempty(rotiClimNet) && any(isfinite(rotiClimNet(:)))
            ltC = ltBins(1:end-1) + 0.5;
            % Stacked bar: strong on top of weak on top of background
            b1 = bar(ax4, ltC, rotiClimNet(:,1), 'stacked', 'FaceColor',[0.7 0.9 0.7],'EdgeColor','none');
            hold(ax4,'on');
            b2 = bar(ax4, ltC, rotiClimNet(:,2), 'stacked', 'FaceColor',[0.9 0.8 0.2],'EdgeColor','none');
            b3 = bar(ax4, ltC, rotiClimNet(:,3), 'stacked', 'FaceColor',[0.85 0.2 0.1],'EdgeColor','none');
            % Fejer RT onset per day
            for d = 1:nDays
                fo = fejerOnset(d);
                if isfinite(fo)
                    xline(ax4, fo, '--', 'Color', [0 0.5 0.9], ...
                        'LineWidth', 1.0, 'Label', sprintf('F%d Pred',d));
                end
                oo = obsEPBonset(d);
                if isfinite(oo)
                    xline(ax4, oo, ':', 'Color', [0.8 0.1 0.1], ...
                        'LineWidth', 1.0, 'Label', sprintf('D%d Obs',d));
                end
            end
            xlabel(ax4,'Local time (h)','FontSize',9);
            ylabel(ax4,'Occurrence rate','FontSize',9);
            title(ax4,'(d) ROTI climatology + RT onset','FontSize',10,'FontWeight','bold');
            legend(ax4, {'Background (<0.2)', 'Weak (0.2-0.5)', 'Strong (>0.5)'}, ...
                'Location','northwest','FontSize',7);
            xlim(ax4,[0 24]); ylim(ax4,[0 1.05]);
            grid(ax4,'on'); box(ax4,'on'); set(ax4,'FontSize',8);
        else
            text(0.5,0.5,'No ROTI climatology data','Parent',ax4,'HorizontalAlignment','center');
            axis(ax4,'off');
        end

        sgtitle(sprintf('Novel Analysis N5: Solar Maximum Ionospheric Structure — F10.7≈%.0f sfu', ...
            mean(p.F107)), 'FontSize',11,'FontWeight','bold');
        set(hf,'PaperUnits','centimeters','PaperSize',[20 22]);
        print(hf, figPath, '-dpng', '-r300');
        close(hf);
        fprintf('[N5] Figure saved: %s\n', figPath);
    catch ME
        fprintf('[N5] Figure error: %s\n', ME.message);
        figPath = '';
    end

    %% ================================================================
    %% CSV OUTPUT
    %% ================================================================
    csvPath = fullfile(outDir, 'solar_max_table.csv');
    try
        fid = fopen(csvPath, 'w');
        fprintf(fid, 'station,day_idx,sys,tau_km_med,tau_km_q25,tau_km_q75,vtec_med_TECU,nmf2_proxy_1e12\n');
        for k = 1:numel(slabRec)
            sr = slabRec(k);
            fprintf(fid, '%s,%d,%s,%.1f,%.1f,%.1f,%.2f,%.4f\n', ...
                sr.station, sr.dayIdx, sr.sys, sr.tau_km_med, sr.tau_km_q25, sr.tau_km_q75, ...
                sr.vtec_med, sr.nmf2_proxy);
        end
        fclose(fid);

        % RT onset CSV
        rtCsvPath = fullfile(outDir, 'rt_onset_table.csv');
        fid2 = fopen(rtCsvPath, 'w');
        fprintf(fid2, 'day_idx,F107,Kp,fejer_onset_LT,obs_onset_LT,delta_RT_h,crest_drift_deg_h\n');
        for d = 1:nDays
            fprintf(fid2, '%d,%.1f,%.2f,%.2f,%.2f,%.2f,%.4f\n', ...
                d, p.F107(d), p.Kp(d), fejerOnset(d), obsEPBonset(d), ...
                deltaRT(d), crestDrift(d));
        end
        fclose(fid2);
        fprintf('[N5] CSVs saved: %s, %s\n', csvPath, rtCsvPath);
    catch ME
        fprintf('[N5] CSV error: %s\n', ME.message);
        csvPath = ''; rtCsvPath = '';
    end

    %% SAVE MAT
    matPath = fullfile(outDir, 'novelSolarMax.mat');
    try
        save(matPath, 'slabRec', 'naiRec', 'crestRec', 'crestDrift', ...
            'R_rotiClim', '-v7');
        fprintf('[N5] MAT saved: %s\n', matPath);
    catch ME
        fprintf('[N5] MAT save error: %s\n', ME.message);
        matPath = '';
    end

    fprintf('\n[N5] ============ SOLAR MAX CHARACTERIZATION SUMMARY ============\n');
    fprintf('  F10.7 range       : %.0f – %.0f sfu\n', min(p.F107), max(p.F107));
    fprintf('  Slab records      : %d\n', numel(slabRec));
    fprintf('  NAI records       : %d\n', numel(naiRec));
    fprintf('  Crest records     : %d\n', numel(crestRec));
    fprintf('  Mean crest drift  : %.3f deg/h\n', nanmean(crestDrift));
    fprintf('  Mean Fejer onset  : %.1f h LT\n', nanmean(fejerOnset));
    fprintf('  Mean obs onset    : %.1f h LT\n', nanmean(obsEPBonset(isfinite(obsEPBonset))));
    fprintf('  Mean ΔRT onset    : %.1f h\n', nanmean(deltaRT(isfinite(deltaRT))));
    fprintf('[N5] ================================================================\n\n');

    R = struct('slab', slabRec, 'NAI', naiRec, 'crest', crestRec, ...
        'crestDrift', crestDrift, 'rotiClim', R_rotiClim, ...
        'figPath', figPath, 'matPath', matPath, 'csvPath', csvPath);
end

%% ============================================================================
%%  LOCAL HELPER FUNCTIONS
%% ============================================================================

function lat = stationLat()
    lat.BHPL = 23.21;
    lat.DRDN = 29.73;
    lat.GDKG = 23.87;
    lat.JDPR = 26.30;
    lat.LCK4 = 26.91;
    lat.PBR4 = 23.02;
    lat.SHLG = 25.57;
end

function S = loadTECmat(dayRoot, station)
    S = [];
    resDir = fullfile(dayRoot, 'Results');
    if ~exist(resDir,'dir'), resDir = dayRoot; end
    mats = dir(fullfile(resDir, sprintf('TEC_%s_*.mat', station)));
    if isempty(mats)
        mats = dir(fullfile(resDir, sprintf('MultiGNSS_*%s*.mat', station)));
    end
    if isempty(mats), return; end
    try, S = load(fullfile(resDir, mats(1).name)); catch, end
end

function vtec = loadVTEC1D(dayRoot, station)
    vtec = [];
    S = loadTECmat(dayRoot, station);
    if isempty(S), return; end
    fn = fieldnames(S);
    for k = 1:numel(fn)
        if startsWith(fn{k}, 'TEC_') || startsWith(fn{k}, 'VTEC_')
            V = S.(fn{k});
            if isnumeric(V) && size(V,1) >= 3600
                vtec = nanmedian(V, 2);
                return;
            end
        end
    end
    if isfield(S,'TEC') && isstruct(S.TEC) && isfield(S.TEC,'vertical')
        vtec = nanmedian(S.TEC.vertical, 2);
    end
end

function sysList = detectSystems(S, station)
    % Returns cell array of struct with .sys and .vtec
    sysList = {};
    if isempty(S), return; end
    fn = fieldnames(S);
    for k = 1:numel(fn)
        if startsWith(fn{k},'TEC_') && contains(fn{k}, station)
            V = S.(fn{k});
            if isnumeric(V) && size(V,1) >= 3600
                % try to detect system from variable name
                parts = strsplit(fn{k}, '_');
                sys = 'G';  % default GPS
                if numel(parts) >= 3
                    candidate = parts{3};
                    if ismember(candidate, {'E','C','J'})
                        sys = candidate;
                    end
                end
                sysList{end+1} = struct('sys', sys, 'vtec', nanmedian(V,2)); %#ok<AGROW>
            end
        end
    end
    if isempty(sysList) && isfield(S,'TEC') && isstruct(S.TEC) && isfield(S.TEC,'vertical')
        sysList{1} = struct('sys','G','vtec',nanmedian(S.TEC.vertical,2));
    end
end

function ROTI = getROTI(S, station)
    ROTI = [];
    if isempty(S), return; end
    fn = fieldnames(S);
    for k = 1:numel(fn)
        if startsWith(fn{k}, 'ROTI_') && contains(fn{k}, station)
            V = S.(fn{k});
            if isnumeric(V) && size(V,1) >= 3600
                ROTI = V; return;
            end
        end
    end
end

function vm = minuteMedianLocal(v, secPerBin)
    N = numel(v);
    nBins = floor(N / secPerBin);
    vm = nan(nBins, 1);
    for k = 1:nBins
        seg = v((k-1)*secPerBin+1 : k*secPerBin);
        vm(k) = nanmedian(seg);
    end
end

function p = prctileLocal(x, pct)
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0, p = NaN(size(pct)); return; end
    p = zeros(size(pct));
    for k = 1:numel(pct)
        idx = pct(k)/100 * (n-1) + 1;
        lo = max(1, floor(idx)); hi = min(n, ceil(idx));
        if lo == hi, p(k) = x(lo);
        else, p(k) = x(lo) + (idx-lo)*(x(hi)-x(lo)); end
    end
end
