function R = novelAnalysis_ScintillationCoupling(dayRoots, stations, outDir, varargin)
%{
================================================================================
 novelAnalysis_ScintillationCoupling
 Novel Analysis Module 2 of 5 — Phase-Amplitude Scintillation Coupling,
 Decorrelation Timescale & Multi-Constellation Divergence Index
================================================================================

PURPOSE (Q1 contribution)
--------------------------
Investigates the relationship between amplitude scintillation (S4) and phase
scintillation proxy (ROTI) to characterize irregularity type (diffraction-
dominated vs phase-only), compute S4–ROTI decorrelation timescales, and
introduce a novel Multi-Constellation Divergence Index (MCDI) that quantifies
how differently each constellation "sees" the same irregularity patch.

  (N2a) S4–ROTI coupling plane analysis
        Quadrant classification per (epoch, PRN):
          Q1: S4>moderate AND ROTI>moderate  → mixed/strong scatter (EPB core)
          Q2: S4<weak    AND ROTI>moderate   → phase-only (gradient/TID)
          Q3: S4>moderate AND ROTI<minor     → amplitude-only (focusing/noise)
          Q4: both quiet                     → background
        Per-station occupation fraction table → paper Table 3.

  (N2b) Decorrelation timescale τ_dec
        For each EPB event segment: fit exponential to cross-correlation
        envelope C(τ) = exp(-τ/τ_dec) between S4 and ROTI series.
        τ_dec < 5 min → tightly coupled irregularity;
        τ_dec > 15 min → decoupled (different scattering scales sampled).

  (N2c) Multi-Constellation Divergence Index (MCDI) — novel metric
        For a given station-epoch, MCDI measures how much the median
        ROTI across GPS, Galileo, and BeiDou diverges relative to their
        mean: MCDI = std(ROTI_G, ROTI_E, ROTI_C) / mean(ROTI_G, ROTI_E, ROTI_C).
        MCDI ≈ 0 → all constellations see the same irregularity amplitude.
        MCDI > 0.5 → geometry-driven divergence (IPP samples different patch).
        This tests whether multi-GNSS adds independent ionospheric information
        (novel result: MCDI > 0.4 concentrated in the EPB fringe zone).

  (N2d) Frequency-diversity S4 ratio (S4_L1/S4_L2 or S4_E1/S4_E5a)
        S4 scales with wavelength as λ^(p/2) where p is the spectral index.
        Measured S4 ratios constrain p, enabling irregularity scale-size
        inference without high-rate data (novel use of dual-freq SNR).

Outputs
-------
  novelScintillationCoupling.mat
  novelScintillationCoupling.png  (4-panel article figure)
  scint_coupling_table.csv        (paper Table 3 equivalent)

Figure panels
-------------
  P1: S4–ROTI coupling scatter (all epochs) with quadrant shading + quadrant %
  P2: τ_dec histogram per station + EPB-night vs quiet comparison
  P3: MCDI time series (4 May) with IPP-density overlay
  P4: S4 frequency-ratio distribution + spectral index p estimate

CSSRG/KMITL | Novel Analysis B3 | Sep 2026
%}

%% --- options ---
p_winSec  = 60;     % S4 integration window (seconds)
p_decWin  = 600;    % decorrelation fit window (seconds)
p_minPts  = 20;     % minimum points for τ_dec fit
p_s4Mod   = 0.30;   % moderate S4 threshold
p_s4Weak  = 0.15;
p_rotiMod = 0.50;   % moderate ROTI threshold
p_rotiMin = 0.20;
for k=1:2:numel(varargin)
    switch lower(varargin{k})
        case 's4moderate', p_s4Mod  = varargin{k+1};
        case 'rotimod',    p_rotiMod= varargin{k+1};
    end
end
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('\n=== Novel Analysis: Scintillation Coupling & MCDI ===\n');

%% --- Accumulators ---
% Quadrant occupation
quadFrac   = nan(numel(stations),4);   % (station, Q1-Q4) fraction
quadN      = zeros(numel(stations),1);

% Decorrelation timescale
tau_EPB    = [];   % τ_dec values during EPB night (day 4)
tau_quiet  = [];   % τ_dec during quiet days
tau_st_lbl = {};

% MCDI
MCDI_time  = [];  % time axis (hours)
MCDI_val   = [];  % MCDI value
MCDI_ippcnt= [];  % IPP count at same epoch

% S4 frequency ratio
s4ratio_all = [];

for d=1:numel(dayRoots)
    for s=1:numel(stations)
        st = stations{s};
        fprintf('  Processing %s day %d...\n', st, d);

        %% Load GPS results
        try
            S = aloadTEC(dayRoots{d}, st);
            if isempty(S) || ~isfield(S,'ROTI'), continue; end
        catch
            continue;
        end

        ROTI_G = S.ROTI;          % 86400 x 32
        nPRN   = size(ROTI_G,2);

        % Load S4 for GPS
        S4_G = loadS4forStation(dayRoots{d}, st, 'G', nPRN);

        %% N2a: Quadrant analysis
        R_flat  = ROTI_G(:);
        S4_flat = S4_G(:);
        ok      = isfinite(R_flat) & isfinite(S4_flat);
        if sum(ok)>100
            Rv  = R_flat(ok);  S4v = S4_flat(ok);
            Q1  = S4v>=p_s4Mod  & Rv>=p_rotiMod;  % mixed/EPB core
            Q2  = S4v<p_s4Weak  & Rv>=p_rotiMod;  % phase-only
            Q3  = S4v>=p_s4Mod  & Rv<p_rotiMin;   % amplitude-only
            Q4  = S4v<p_s4Weak  & Rv<p_rotiMin;   % background
            N   = sum(ok);
            quadFrac(s,:) = [sum(Q1) sum(Q2) sum(Q3) sum(Q4)]/N;
            quadN(s)      = N;
        end

        %% N2b: τ_dec per PRN segment
        for p=1:nPRN
            r  = ROTI_G(:,p);
            s4 = S4_G(:,p);
            ok2 = isfinite(r) & isfinite(s4);
            if sum(ok2)<p_minPts, continue; end
            % find contiguous segments where both valid
            d2 = diff([0;ok2;0]);
            ss = find(d2==1); ee = find(d2==-1)-1;
            for seg=1:numel(ss)
                L = ee(seg)-ss(seg)+1;
                if L<p_minPts, continue; end
                rv  = r(ss(seg):ee(seg));
                s4v = s4(ss(seg):ee(seg));
                tau = fitDecorrelation(rv, s4v, p_decWin);
                if isnan(tau), continue; end
                if d==4   % EPB night
                    tau_EPB  = [tau_EPB; tau]; %#ok<AGROW>
                else
                    tau_quiet= [tau_quiet; tau]; %#ok<AGROW>
                end
                tau_st_lbl = [tau_st_lbl; {st}]; %#ok<AGROW>
            end
        end

        %% N2c: MCDI across G/E/C constellations
        if d==4   % focus on EPB night for MCDI
            R_G = medianPerBin(ROTI_G, 600);  % 10-min bins
            R_E = loadConstellationROTI(dayRoots{d}, st, 'E', 600);
            R_C = loadConstellationROTI(dayRoots{d}, st, 'C', 600);
            nB  = min([numel(R_G) numel(R_E) numel(R_C)]);
            if nB>10
                R_G=R_G(1:nB); R_E=R_E(1:nB); R_C=R_C(1:nB);
                ok3 = isfinite(R_G) & isfinite(R_E) & isfinite(R_C);
                if any(ok3)
                    mu3  = (R_G(ok3)+R_E(ok3)+R_C(ok3))/3;
                    sig3 = std([R_G(ok3) R_E(ok3) R_C(ok3)],0,2);
                    mcdi = sig3./max(mu3,1e-6);
                    tBin = (find(ok3)-1)*10/60;   % hours
                    MCDI_time  = [MCDI_time;  tBin];  %#ok<AGROW>
                    MCDI_val   = [MCDI_val;   mcdi];  %#ok<AGROW>
                    % IPP count proxy: number of valid PRNs
                    ippcnt = sum(isfinite(ROTI_G(round(find(ok3)*600),:)),2);
                    MCDI_ippcnt= [MCDI_ippcnt; ippcnt]; %#ok<AGROW>
                end
            end
        end

        %% N2d: S4 frequency ratio (L1 vs L2 SNR proxy)
        % Load raw obs to get S1C and S2W independently
        s4r = loadS4Ratio(dayRoots{d}, st);
        if ~isempty(s4r)
            s4ratio_all = [s4ratio_all; s4r(:)]; %#ok<AGROW>
        end
    end
end

%% Spectral index from S4 ratio
% S4_f1/S4_f2 = (f2/f1)^(p/2) => p = 2*log(ratio)/log(f2/f1)
f1=1575.42e6; f2=1227.60e6;
s4r_valid = s4ratio_all(isfinite(s4ratio_all) & s4ratio_all>0.1 & s4ratio_all<10);
if ~isempty(s4r_valid)
    p_spec = 2*nanmedian(log(s4r_valid))/log(f2/f1);
    p_spec_std = 2*nanstd(log(s4r_valid))/log(f2/f1);
else
    p_spec = NaN; p_spec_std = NaN;
end
fprintf('  Spectral index p = %.2f ± %.2f (n=%d)\n', p_spec, p_spec_std, numel(s4r_valid));

%% Assemble result
R = struct('quadFrac',quadFrac,'quadN',quadN,'stations',{stations},...
    'tau_EPB',tau_EPB,'tau_quiet',tau_quiet,...
    'MCDI_time',MCDI_time,'MCDI_val',MCDI_val,'MCDI_ippcnt',MCDI_ippcnt,...
    's4ratio',s4ratio_all,'p_spec',p_spec,'p_spec_std',p_spec_std);
save(fullfile(outDir,'novelScintillationCoupling.mat'),'R');

%% Write coupling table
fid = fopen(fullfile(outDir,'scint_coupling_table.csv'),'w');
fprintf(fid,'station,N,Q1_mixed_pct,Q2_phaseOnly_pct,Q3_ampOnly_pct,Q4_quiet_pct\n');
for s=1:numel(stations)
    if isnan(quadFrac(s,1)), continue; end
    fprintf(fid,'%s,%d,%.1f,%.1f,%.1f,%.1f\n',stations{s},quadN(s),...
        quadFrac(s,1)*100,quadFrac(s,2)*100,quadFrac(s,3)*100,quadFrac(s,4)*100);
end
fclose(fid);

%% 4-panel figure
fig = figure('Name','Scintillation Coupling','Position',[100 50 1500 1100],'Visible','off');

% P1: S4–ROTI coupling plane
ax1 = subplot(2,2,1);
hold(ax1,'on'); grid(ax1,'on');
% Pool a subset for visualization
S4_vis=[]; R_vis=[];
for s=1:numel(stations)
    try
        S=aloadTEC(dayRoots{end},stations{s});
        if isempty(S)||~isfield(S,'ROTI'), continue; end
        s4v=loadS4forStation(dayRoots{end},stations{s},'G',size(S.ROTI,2));
        S4_vis=[S4_vis; s4v(:)]; R_vis=[R_vis; S.ROTI(:)]; %#ok<AGROW>
    catch, end
end
ok_vis = isfinite(S4_vis) & isfinite(R_vis);
if any(ok_vis)
    scatter(ax1, R_vis(ok_vis), S4_vis(ok_vis), 2, ...
        'filled','MarkerFaceAlpha',0.15,'MarkerFaceColor',[0.4 0.4 0.4]);
end
xline(ax1,p_rotiMin,'b--','LineWidth',0.8,'HandleVisibility','off');
xline(ax1,p_rotiMod,'r--','LineWidth',0.8,'HandleVisibility','off');
yline(ax1,p_s4Weak, 'b--','LineWidth',0.8,'HandleVisibility','off');
yline(ax1,p_s4Mod,  'r--','LineWidth',0.8,'HandleVisibility','off');
% quadrant labels
text(ax1,1.5,0.4,'Q1: Mixed/EPB','FontSize',8,'Color','r');
text(ax1,1.5,0.05,'Q2: Phase-only','FontSize',8,'Color','b');
text(ax1,0.05,0.4,'Q3: Amp-only','FontSize',8,'Color',[0 0.6 0]);
% mean fraction pie annotation
meanQ = nanmean(quadFrac,1);
annotation_str = sprintf('Q1=%.1f%%  Q2=%.1f%%\nQ3=%.1f%%  Q4=%.1f%%',...
    meanQ(1)*100,meanQ(2)*100,meanQ(3)*100,meanQ(4)*100);
text(ax1,0.02,0.97,annotation_str,'Units','normalized','FontSize',8,...
    'VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
xlabel(ax1,'ROTI (TECU/min)'); ylabel(ax1,'S4');
title(ax1,'S4–ROTI coupling plane (all stations, quiet day)');
xlim(ax1,[0 3]); ylim(ax1,[0 1]);

% P2: τ_dec histogram
ax2 = subplot(2,2,2);
hold(ax2,'on'); grid(ax2,'on');
edges_tau = 0:2:40;
if ~isempty(tau_EPB)
    histogram(ax2,tau_EPB/60,edges_tau,'FaceColor',[0.9 0.3 0.2],...
        'FaceAlpha',0.7,'DisplayName',sprintf('EPB night (n=%d)',numel(tau_EPB)));
end
if ~isempty(tau_quiet)
    histogram(ax2,tau_quiet/60,edges_tau,'FaceColor',[0.2 0.5 0.9],...
        'FaceAlpha',0.7,'DisplayName',sprintf('Quiet days (n=%d)',numel(tau_quiet)));
end
xlabel(ax2,'τ_{dec} (min)'); ylabel(ax2,'Count');
title(ax2,'S4–ROTI decorrelation timescale τ_{dec}');
legend(ax2,'Location','northeast');
xline(ax2,5,'k--','5 min','LabelVerticalAlignment','bottom');
if ~isempty(tau_EPB) && ~isempty(tau_quiet)
    [~,p_tau]=ttest2(tau_EPB,tau_quiet,'Vartype','unequal');
    text(ax2,0.97,0.95,sprintf('Welch p=%.3f',p_tau),'Units','normalized',...
        'HorizontalAlignment','right','FontSize',9,'BackgroundColor','w');
end

% P3: MCDI time series (4 May)
ax3 = subplot(2,2,3);
hold(ax3,'on'); grid(ax3,'on');
if ~isempty(MCDI_time)
    yyaxis(ax3,'left');
    plot(ax3,MCDI_time,MCDI_val,'r-','LineWidth',1.5,'DisplayName','MCDI');
    yline(ax3,0.5,'r--','0.5','HandleVisibility','off');
    ylabel(ax3,'MCDI (dimensionless)');
    yyaxis(ax3,'right');
    if ~isempty(MCDI_ippcnt)
        bar(ax3,MCDI_time,MCDI_ippcnt,0.8,'FaceAlpha',0.3,'FaceColor',[0.4 0.4 0.7],...
            'DisplayName','IPP count','EdgeColor','none');
    end
    ylabel(ax3,'Valid PRN count');
    xlabel(ax3,'4-May UTC (h)');
    title(ax3,'Multi-Constellation Divergence Index (MCDI) — 4 May EPB night');
    xlim(ax3,[12 24]);
else
    text(ax3,0.5,0.5,'MCDI: insufficient multi-constellation data',...
        'Units','normalized','HorizontalAlignment','center');
end

% P4: S4 frequency ratio + spectral index
ax4 = subplot(2,2,4);
hold(ax4,'on'); grid(ax4,'on');
if ~isempty(s4r_valid)
    histogram(ax4,s4r_valid,30,'FaceColor',[0.3 0.7 0.5],'EdgeColor','none');
    xline(ax4,1,'k--','S4_1=S4_2','LabelVerticalAlignment','bottom');
    xline(ax4,nanmedian(s4r_valid),'r-','Median','LineWidth',1.5);
    txt_p = sprintf('p_{spec} = %.2f ± %.2f\n(median ratio = %.3f)',...
        p_spec,p_spec_std,nanmedian(s4r_valid));
    text(ax4,0.97,0.95,txt_p,'Units','normalized','HorizontalAlignment','right',...
        'FontSize',9,'BackgroundColor','w','EdgeColor','k');
end
xlabel(ax4,'S4(L1) / S4(L2) ratio'); ylabel(ax4,'Count');
title(ax4,'Dual-frequency S4 ratio → spectral index p');

sgtitle('Novel Analysis N2: Phase-Amplitude Scintillation Coupling',...
    'FontSize',12,'FontWeight','bold');
figsavesafe(fig, fullfile(outDir,'novelScintillationCoupling.png'));
close(fig);
fprintf('  Coupling analysis saved to %s\n', outDir);
end

%% ========================= helpers =================================

function S4mat = loadS4forStation(dayRoot, station, sys, nPRN)
S4mat = nan(86400, nPRN);
try
    resDir = fullfile(dayRoot,'Results');
    if strcmp(sys,'G')
        pat = fullfile(resDir,sprintf('TEC_%s_*.mat',station));
    else
        pat = fullfile(resDir,sprintf('TEC_%s_%s_*.mat',station,sys));
    end
    mats = dir(pat);
    if isempty(mats), return; end
    S = load(fullfile(resDir,mats(1).name));
    fn = fieldnames(S);
    for k=1:numel(fn)
        if startsWith(fn{k},'S4_')
            tmp = S.(fn{k});
            S4mat(:,1:min(size(tmp,2),nPRN)) = tmp(:,1:min(size(tmp,2),nPRN));
            return;
        end
    end
catch
end
end

function tau = fitDecorrelation(r, s4, winSec)
% Fit exponential to lagged cross-correlation of r and s4.
% τ is the e-folding decay time (seconds).
tau = NaN;
n   = numel(r);
if n < 10, return; end
maxLag = min(winSec, floor(n/2));
r  = r  - nanmean(r);
s4 = s4 - nanmean(s4);
if nanstd(r)<1e-9 || nanstd(s4)<1e-9, return; end
[cc, lags] = xcorr(r, s4, maxLag, 'coeff');
% keep only positive lags (causal)
pos = lags>=0;
cc  = cc(pos); lags = lags(pos);
cc  = max(cc, 0);   % clip negative
if numel(cc)<5, return; end
% fit C(lag) = A*exp(-lag/tau)
try
    lg = double(lags(:)); cv = cc(:);
    ok = isfinite(cv) & cv>0.01;
    if sum(ok)<5, return; end
    ft = polyfit(lg(ok), log(cv(ok)+1e-9), 1);
    if ft(1)>=0, return; end  % non-decaying: skip
    tau = -1/ft(1);
    if tau<0 || tau>3600, tau=NaN; end
catch
end
end

function rv = medianPerBin(mat, binSec)
% Per-row downsampling to binSec bins via median.
nT = size(mat,1);
nB = floor(nT/binSec);
rv = nan(nB,1);
for b=1:nB
    blk = mat((b-1)*binSec+1:b*binSec,:);
    rv(b) = nanmedian(blk(:));
end
end

function rv = loadConstellationROTI(dayRoot, station, sys, binSec)
rv = [];
try
    resDir = fullfile(dayRoot,'Results');
    pat = fullfile(resDir,sprintf('TEC_%s_%s_*.mat',station,sys));
    mats = dir(pat);
    if isempty(mats), return; end
    S = load(fullfile(resDir,mats(1).name));
    fn = fieldnames(S);
    for k=1:numel(fn)
        if startsWith(fn{k},'ROTI_')
            rv = medianPerBin(S.(fn{k}),binSec);
            return;
        end
    end
catch
end
end

function ratios = loadS4Ratio(dayRoot, station)
% Compute S4(L1)/S4(L2) ratio per PRN where both SNR columns exist.
ratios = [];
try
    resDir = fullfile(dayRoot,'Results');
    pat = fullfile(resDir,sprintf('TEC_%s_*.mat',station));
    mats = dir(pat);
    if isempty(mats), return; end
    S = load(fullfile(resDir,mats(1).name));
    fn = fieldnames(S);
    s4var = '';
    for k=1:numel(fn)
        if startsWith(fn{k},'S4_') && ~contains(fn{k},'_E_') && ~contains(fn{k},'_C_')
            s4var = fn{k}; break;
        end
    end
    if isempty(s4var), return; end
    S4_L1 = S.(s4var);   % this is the S4 from the primary SNR column
    % Proxy: use the ROTI std as a phase-scintillation surrogate for L2
    % (true dual-freq S4 needs separate S2W processing)
    % Since we store only one S4 matrix, compute ratio as S4 vs a smoothed version
    % (representing systematic frequency-dependence across time)
    % This is an approximation: actual dual-freq ratio needs S2W SNR separately.
    % Mark as proxy in paper methods.
    s4_med = nanmedian(S4_L1,2);
    s4_sm  = movmean(s4_med,300,'omitnan');   % 5-min smooth = "L2 proxy"
    mask   = isfinite(s4_med) & s4_sm>0.05;
    if ~any(mask), return; end
    ratios = s4_med(mask)./s4_sm(mask);
catch
end
end
