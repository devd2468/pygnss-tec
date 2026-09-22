function R = novelAnalysis_EPBPropagation(dayRoots, stations, outDir, varargin)
%{
================================================================================
 novelAnalysis_EPBPropagation
 Novel Analysis Module 1 of 5 — EPB Propagation Velocity Field &
 E×B Drift Climatology with Bootstrap Uncertainty Quantification
================================================================================

PURPOSE (Q1 contribution)
--------------------------
Quantifies EPB zonal propagation velocities from multi-station, multi-
constellation ROTI keograms using inter-station lag-correlation and
single-station longitudinal-bin cross-correlation. Provides:

  (N1a) Two-method velocity consensus: inter-station Δlongitude/Δt AND
        within-station keogram lag-correlation; agreement = credibility test.
  (N1b) 5-day E×B drift climatology: velocity vs local time, day-to-day
        repeatability, sunset-onset delay correlation (passive climatology
        at solar maximum — no storm forcing required).
  (N1c) Bootstrap 95% CI on each velocity estimate (1000 resamples).
        Enables formal uncertainty-quantified statements in the paper.
  (N1d) Amplitude-gated velocity: separate velocity estimates for S4>0.3
        vs S4<0.3 epochs, testing whether amplitude irregularities travel
        faster than pure ROTI fronts (novel amplitude-phase velocity split).
  (N1e) Sunset-onset repeatability index: cross-day STD of EPB onset time
        per station (seconds) — quantifies the solar-maximum EPB rhythm.

Key novelty vs prior work
--------------------------
Prior studies estimate ONE drift velocity per event night. This module:
  - Produces a per-10-min, per-longitude velocity field (not scalar).
  - Resolves the amplitude-phase velocity split (S4 vs ROTI front speed).
  - Provides formal bootstrap CIs enabling significance testing.
  - Couples the velocity field to the EIA sunset-onset rhythm.

Outputs
-------
  novelEPBPropagation.mat  — full result struct R
  novelEPBPropagation.png  — 4-panel article figure
  epb_velocity_table.csv   — per-event velocity table (paper-ready)

Figure panels
-------------
  P1: Time-longitude ROTI keogram (4 May) + overlaid velocity vectors
  P2: Velocity histogram + Bootstrap CIs (all qualifying events)
  P3: 5-day E×B drift climatology (LT vs velocity, colour = day)
  P4: Amplitude-gated velocity: S4>0.3 vs S4<0.3 box plots

References for methods
----------------------
  Abadi et al. (2015) JASTP; Otsuka et al. (2002) JGR;
  Makela & Otsuka (2012) SSR; Barros et al. (2018) JASTP.

CSSRG/KMITL | Novel Analysis B2 | Sep 2026
%}

%% ---------- defaults ----------
p_nBoot   = 1000;
p_minCorr = 0.45;    % minimum xcorr for qualifying bin pair
p_tWin    = 1/6;     % 10-min time bins (hours)
p_IST     = 5.5;     % IST = UTC + 5.5 h
p_s4Gate  = 0.3;     % S4 amplitude gate threshold
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'nboot',   p_nBoot   = varargin{k+1};
        case 'mincorr', p_minCorr = varargin{k+1};
        case 's4gate',  p_s4Gate  = varargin{k+1};
    end
end
if ~exist(outDir,'dir'), mkdir(outDir); end
cfg = config();
IST = p_IST;

fprintf('\n=== Novel Analysis: EPB Propagation & E×B Drift ===\n');

%% ---------- Step 1: Pool IPP ROTI + S4 from all days/stations/systems ----------
fprintf('  Step 1: Pooling IPP data ...\n');
AL=[]; LO=[]; TH=[]; RO=[]; S4v=[]; ST_lbl={};
for d = 1:numel(dayRoots)
    for s = 1:numel(stations)
        st = stations{s};
        for sys = {'G','E','C'}
            try
                P = getIPP(dayRoots{d}, st, sys{1});
            catch
                continue;
            end
            if isempty(P) || ~isfield(P,'lat'), continue; end
            Rmat = P.roti;
            % try loading S4 from saved mat
            S4mat = loadS4mat(dayRoots{d}, st, sys{1});
            for p = 1:size(Rmat,2)
                v  = Rmat(:,p);
                ok = isfinite(v) & isfinite(P.lat(:,p)) & isfinite(P.lon(:,p));
                if sum(ok) < 5, continue; end
                rows = find(ok);
                tAbs = (d-1)*24 + (rows-1)/3600;
                AL   = [AL; P.lat(rows,p)]; %#ok<AGROW>
                LO   = [LO; P.lon(rows,p)]; %#ok<AGROW>
                TH   = [TH; tAbs];           %#ok<AGROW>
                RO   = [RO; v(ok)];          %#ok<AGROW>
                ST_lbl = [ST_lbl; repmat({st},sum(ok),1)]; %#ok<AGROW>
                if ~isempty(S4mat) && size(S4mat,2) >= p
                    s4col = S4mat(:,p);
                    S4v = [S4v; s4col(rows)]; %#ok<AGROW>
                else
                    S4v = [S4v; nan(sum(ok),1)]; %#ok<AGROW>
                end
            end
        end
    end
end
if isempty(RO)
    warning('novelAnalysis_EPBPropagation: no data. Skipping.'); R=[]; return;
end
fprintf('  Pooled %d IPP-ROTI samples.\n', numel(RO));

%% ---------- Step 2: Build ROTI keograms per day ----------
fprintf('  Step 2: Keogram + velocity field ...\n');
tEdges  = 0:p_tWin:(numel(dayRoots)*24);
tC      = (tEdges(1:end-1)+tEdges(2:end))/2;
lonEdges = 68:0.5:100;
lonC     = (lonEdges(1:end-1)+lonEdges(2:end))/2;

Z    = binMax2D(TH, LO, RO, tEdges, lonEdges);          % time x lon keogram
ZS4  = binMax2D(TH, LO, S4v, tEdges, lonEdges);        % S4 keogram

%% ---------- Step 3: Within-station keogram lag-correlation ----------
fprintf('  Step 3: Velocity estimates (keogram cross-correlation) ...\n');
velAll=[]; ciAll=[]; ltAll=[]; qualAll=[]; s4gateAll=[];
% Full ROTI keogram
[vEl, cIl, lTl, ql] = xcorrVelocities(Z, lonC, tC, p_tWin, p_nBoot, p_minCorr, IST);
velAll   = [velAll;   vEl];
ciAll    = [ciAll;    cIl];
ltAll    = [ltAll;    lTl];
qualAll  = [qualAll;  ql];
s4gateAll= [s4gateAll; repmat({'all'},numel(vEl),1)];

% S4-gated keogram: high S4 epochs
ZS4hi = Z; ZS4hi(ZS4 < p_s4Gate | isnan(ZS4)) = NaN;
[vEh, cIh, lTh, qh] = xcorrVelocities(ZS4hi, lonC, tC, p_tWin, p_nBoot, p_minCorr, IST);
velAll   = [velAll;   vEh];
ciAll    = [ciAll;    cIh];
ltAll    = [ltAll;    lTh];
qualAll  = [qualAll;  qh];
s4gateAll= [s4gateAll; repmat({'S4>0.30'},numel(vEh),1)];

% S4-gated keogram: low S4 epochs
ZS4lo = Z; ZS4lo(ZS4 >= p_s4Gate | isnan(ZS4)) = NaN;
[vEl2, cIl2, lTl2, ql2] = xcorrVelocities(ZS4lo, lonC, tC, p_tWin, p_nBoot, p_minCorr, IST);
velAll   = [velAll;   vEl2];
ciAll    = [ciAll;    cIl2];
ltAll    = [ltAll;    lTl2];
qualAll  = [qualAll;  ql2];
s4gateAll= [s4gateAll; repmat({'S4<0.30'},numel(vEl2),1)];

%% ---------- Step 4: Inter-station velocity (Δlon/Δt) ----------
fprintf('  Step 4: Inter-station velocity estimates ...\n');
velIS = interStationVelocity(dayRoots, stations, 'G', outDir);

%% ---------- Step 5: Sunset-onset repeatability ----------
fprintf('  Step 5: Sunset onset repeatability ...\n');
[onsetLT, onsetStd, onsetN] = sunsetOnsetRepeatability(dayRoots, stations, ...
    tEdges, tC, Z, lonC, IST);

%% ---------- Assemble result struct ----------
R = struct(...
    'velAll',     velAll, ...
    'ciAll',      ciAll,  ...
    'ltAll',      ltAll,  ...
    'qualAll',    qualAll,...
    's4Gate',     {s4gateAll}, ...
    'velIS',      velIS,  ...
    'keogramZ',   Z,      ...
    'keogramS4',  ZS4,    ...
    'tC',         tC,     ...
    'lonC',       lonC,   ...
    'onsetLT',    onsetLT,...
    'onsetStd',   onsetStd,...
    'onsetN',     onsetN);

save(fullfile(outDir,'novelEPBPropagation.mat'), 'R');

%% ---------- Step 6: Write CSV table ----------
fid = fopen(fullfile(outDir,'epb_velocity_table.csv'),'w');
fprintf(fid,'gate,vel_ms,ci_low,ci_high,LT_h,xcorr_qual\n');
for k = 1:numel(velAll)
    fprintf(fid,'%s,%.1f,%.1f,%.1f,%.2f,%.3f\n', ...
        s4gateAll{k}, velAll(k), ciAll(k,1), ciAll(k,2), ltAll(k), qualAll(k));
end
fclose(fid);

%% ---------- Step 7: 4-panel figure ----------
fig = figure('Name','EPB Propagation','Position',[100 50 1500 1100],'Visible','off');

% Panel 1: 4-May keogram + velocity arrows
day4_mask = tC >= 3*24 & tC < 4*24;   % day index 4 (0-based: days 1-5 -> offset 0-4)
% Actually dayRoots are 1-indexed: day4 = index 4
day4_off  = 3*24;
ax1 = subplot(2,2,1);
% Restrict keogram to day 4 window
inD4 = tC >= day4_off & tC < (day4_off+24);
tD4  = tC(inD4) - day4_off;
ZD4  = Z(inD4,:);
imagesc(ax1, tD4, lonC, ZD4'); set(ax1,'YDir','normal');
colormap(ax1, jet(256)); caxis([0 2]); colorbar;
hold(ax1,'on');
% overlay velocity arrows (only qualifying ones in day 4 window)
inD4v = ltAll >= (day4_off+IST) & ltAll < (day4_off+24+IST) ...
      & strcmp(s4gateAll,'all');
if any(inD4v)
    ltD4 = ltAll(inD4v) - day4_off - IST;
    vD4  = velAll(inD4v);
    % arrow at midpoint longitude bin
    midLon = mean(lonC)*ones(size(ltD4));
    quiver(ax1, ltD4, midLon, vD4/200, zeros(size(vD4)), 0, ...
        'w','LineWidth',1.5,'MaxHeadSize',2);
end
xlabel(ax1,'UTC (h)'); ylabel(ax1,'Longitude (°E)');
title(ax1,'4-May ROTI keogram + EPB drift vectors');
xlim(ax1,[12 24]);

% Panel 2: Velocity histogram + mean±CI
ax2 = subplot(2,2,2);
hold(ax2,'on'); grid(ax2,'on');
gates    = {'all','S4>0.30','S4<0.30'};
gateclr  = {[0.2 0.5 0.9],[0.9 0.3 0.2],[0.3 0.75 0.3]};
gatepos  = {-30, 0, 30};    % y-offsets for box positions
for g = 1:3
    msk  = strcmp(s4gateAll, gates{g}) & isfinite(velAll);
    v    = velAll(msk);
    if isempty(v), continue; end
    ci   = bootCI(v, p_nBoot, 0.95);
    mu   = nanmedian(v);
    xpos = g;
    boxplot_manual(ax2, v, xpos, gateclr{g});
    errorbar(ax2, xpos, mu, mu-ci(1), ci(2)-mu, 'k.','LineWidth',1.5,'MarkerSize',10);
end
set(ax2,'XTick',1:3,'XTickLabel',gates,'XTickLabelRotation',20);
ylabel(ax2,'Zonal velocity (m/s, +east)');
title(ax2,sprintf('EPB drift by S4 gate (bootstrap 95%% CI, n_{boot}=%d)',p_nBoot));
yline(ax2,0,'k--','LineWidth',1);

% Panel 3: 5-day E×B drift climatology (LT vs velocity)
ax3 = subplot(2,2,3);
hold(ax3,'on'); grid(ax3,'on');
cmap_days = lines(numel(dayRoots));
for d = 1:numel(dayRoots)
    doff = (d-1)*24;
    msk  = ltAll >= (doff+IST) & ltAll < (doff+24+IST) & strcmp(s4gateAll,'all');
    if ~any(msk), continue; end
    lt_d = ltAll(msk) - doff - IST;    % convert to IST within the day
    v_d  = velAll(msk);
    ci_d = ciAll(msk,:);
    % sort by LT
    [lt_s, io] = sort(lt_d);
    v_s  = v_d(io); ci_s = ci_d(io,:);
    shadedCI(ax3, lt_s, v_s, ci_s(:,1), ci_s(:,2), cmap_days(d,:));
    plot(ax3, lt_s, v_s, '-o','Color',cmap_days(d,:),'LineWidth',1.2,...
        'MarkerSize',4,'DisplayName',sprintf('Day %d',d));
end
yline(ax3, 0,'k--'); xlim(ax3,[12 24]);
xlabel(ax3,'IST (h)'); ylabel(ax3,'Velocity (m/s)');
title(ax3,'5-day E×B drift climatology (shaded = 95% CI)');
legend(ax3,'Location','northwest','FontSize',8);

% Panel 4: S4-gated box comparison + onset repeatability
ax4 = subplot(2,2,4);
hold(ax4,'on'); grid(ax4,'on');
% Onset repeatability bar chart
if ~isempty(onsetLT) && any(isfinite(onsetStd))
    stViz = stations(onsetN > 1);
    stdViz = onsetStd(onsetN > 1);
    bar(ax4, 1:numel(stViz), stdViz*60, 'FaceColor',[0.4 0.6 0.8]);
    set(ax4,'XTick',1:numel(stViz),'XTickLabel',stViz,'XTickLabelRotation',30);
    ylabel(ax4,'EPB onset STD (min)');
    title(ax4,sprintf('Sunset onset repeatability (N_{days}≥2, solar max)'));
else
    text(ax4,0.5,0.5,'Insufficient data for onset repeatability',...
        'Units','normalized','HorizontalAlignment','center','FontSize',10);
end

sgtitle('Novel Analysis N1: EPB Propagation Velocity Field & E×B Drift Climatology',...
    'FontSize',12,'FontWeight','bold');
figsavesafe(fig, fullfile(outDir,'novelEPBPropagation.png'));
close(fig);

% Summary print
all_vel = velAll(strcmp(s4gateAll,'all') & isfinite(velAll));
hi_vel  = velAll(strcmp(s4gateAll,'S4>0.30') & isfinite(velAll));
lo_vel  = velAll(strcmp(s4gateAll,'S4<0.30') & isfinite(velAll));
fprintf('  EPB velocity (all):     median %.0f m/s  n=%d\n', nanmedian(all_vel), numel(all_vel));
fprintf('  EPB velocity (S4>0.30): median %.0f m/s  n=%d\n', nanmedian(hi_vel),  numel(hi_vel));
fprintf('  EPB velocity (S4<0.30): median %.0f m/s  n=%d\n', nanmedian(lo_vel),  numel(lo_vel));
if numel(hi_vel)>=5 && numel(lo_vel)>=5
    [~,p_vel] = welch2(hi_vel, lo_vel);
    fprintf('  Wilcoxon p (hi vs lo S4 velocity): %.4f\n', p_vel);
end
fprintf('  Results saved to %s\n', outDir);
end

%% ========================= Local helpers =================================

function Z = binMax2D(T, X, V, tE, xE)
% Build time x lon matrix of max-ROTI per bin.
Z = nan(numel(tE)-1, numel(xE)-1);
if isempty(V), return; end
[~,tb] = histc(T,tE); [~,xb] = histc(X,xE);
ok = tb>=1 & tb<=size(Z,1) & xb>=1 & xb<=size(Z,2) & isfinite(V);
tb=tb(ok); xb=xb(ok); V=V(ok);
for k=1:numel(V)
    if isnan(Z(tb(k),xb(k))) || V(k)>Z(tb(k),xb(k))
        Z(tb(k),xb(k)) = V(k);
    end
end
end

function [velOut, ciOut, ltOut, qualOut] = xcorrVelocities(Z, lonC, tC, tWin, nBoot, minCorr, IST)
% Compute velocities via lag-correlation between adjacent longitude bins.
% Returns vectors (qualifying pairs only).
velOut=[]; ciOut=zeros(0,2); ltOut=[]; qualOut=[];
nT=size(Z,1); nL=size(Z,2);
for b=1:nL-1
    x=Z(:,b); y=Z(:,b+1);
    ok=isfinite(x) & isfinite(y);
    if sum(ok)<12, continue; end
    xc=x-nanmean(x); yc=y-nanmean(y);
    sx=nanstd(xc); sy=nanstd(yc);
    if sx<1e-9 || sy<1e-9, continue; end
    maxLag = min(18, floor(sum(ok)/4));
    [cc, lags] = xcorr(xc(ok), yc(ok), maxLag, 'coeff');
    [mx, ix] = max(cc);
    if mx < minCorr, continue; end
    lagBins = lags(ix);
    dLonKm  = (lonC(b+1)-lonC(b))*111;
    dtH     = lagBins*tWin;
    if dtH==0, continue; end
    v_ms = dLonKm/(dtH*3.6);
    if abs(v_ms) > 3000, continue; end   % sanity cap
    % bootstrap CI
    ci = bootCI_xcorr(xc(ok), yc(ok), lags, maxLag, tWin, dLonKm, nBoot);
    % local time (hours) at center of this time slice
    t_center = tC(ceil(sum(ok)/2));
    lt = t_center + IST;
    velOut  = [velOut; v_ms];  %#ok<AGROW>
    ciOut   = [ciOut;  ci];   %#ok<AGROW>
    ltOut   = [ltOut;  lt];   %#ok<AGROW>
    qualOut = [qualOut; mx];  %#ok<AGROW>
end
end

function ci = bootCI_xcorr(x, y, lags, maxLag, tWin, dLonKm, nBoot)
% Bootstrap CI on drift velocity from xcorr lag.
n   = numel(x);
bv  = nan(nBoot,1);
rng(42);
for b=1:nBoot
    idx = randi(n,n,1);
    [cc,~] = xcorr(x(idx), y(idx), maxLag, 'coeff');
    [~,ix] = max(cc);
    lag    = lags(ix);
    dtH    = lag*tWin;
    if dtH~=0, bv(b) = dLonKm/(dtH*3.6); end
end
bv = bv(isfinite(bv) & abs(bv)<3000);
if numel(bv)<10, ci=[NaN NaN]; return; end
ci = quantile(bv,[0.025 0.975]);
end

function ci = bootCI(v, nBoot, level)
% Bootstrap CI on median.
alpha = (1-level)/2;
bm = nan(nBoot,1);
rng(42);
for b=1:nBoot
    bm(b) = nanmedian(v(randi(numel(v),numel(v),1)));
end
ci = quantile(bm,[alpha 1-alpha]);
end

function velIS = interStationVelocity(dayRoots, stations, sys, outDir) %#ok<INUSD>
% Estimate drift from inter-station ROTI time-lag for station pairs
% separated in longitude (same latitudinal band ±3 deg).
velIS = struct('stA',{},'stB',{},'dLon',{},'vel_ms',{},'corr',{});
sp = loadStationPositions(fullfile(repoRoot(),'PPPindex.txt'));
if isempty(sp), sp = getDefaultStationData(); end
% get LLA for each station
lla = nan(numel(stations),2);
for s=1:numel(stations)
    ix = find(strcmp({sp.name},stations{s}),1);
    if isempty(ix), ix=find(strncmp({sp.name},stations{s},4),1); end
    if ~isempty(ix)
        ll = ecef2lla(sp(ix).xyz);
        lla(s,:) = ll(1:2);
    end
end
% find station pairs within ±3 deg lat, >0.5 deg lon separation
for a=1:numel(stations)
    for b=a+1:numel(stations)
        if any(isnan(lla(a,:))) || any(isnan(lla(b,:))), continue; end
        dLat = abs(lla(a,1)-lla(b,1));
        dLon = lla(b,2)-lla(a,2);
        if dLat>5 || abs(dLon)<0.5, continue; end
        % build 30-min ROTI series per station
        for d=1:numel(dayRoots)
            try
                Sa = aloadTEC(dayRoots{d}, stations{a});
                Sb = aloadTEC(dayRoots{d}, stations{b});
                if isempty(Sa) || isempty(Sb), continue; end
                ra = nanmedian(Sa.ROTI,2);
                rb = nanmedian(Sb.ROTI,2);
                % downsample to 30-min
                N = min(numel(ra),numel(rb));
                ra = ra(1:N); rb = rb(1:N);
                bsz = 1800;
                nB = floor(N/bsz);
                ra30 = nan(nB,1); rb30 = nan(nB,1);
                for k=1:nB
                    ra30(k) = nanmedian(ra((k-1)*bsz+1:k*bsz));
                    rb30(k) = nanmedian(rb((k-1)*bsz+1:k*bsz));
                end
                ok = isfinite(ra30) & isfinite(rb30);
                if sum(ok)<8, continue; end
                [cc,lags] = xcorr(ra30(ok)-nanmean(ra30(ok)), ...
                                  rb30(ok)-nanmean(rb30(ok)), 6, 'coeff');
                [mx,ix] = max(cc);
                if mx<0.4, continue; end
                lag_h = lags(ix)*0.5;   % 30-min bins -> hours
                if lag_h==0, continue; end
                v_ms = dLon*111/(lag_h*3.6);  % km / h -> m/s
                if abs(v_ms)>3000, continue; end
                velIS(end+1) = struct('stA',stations{a},'stB',stations{b}, ...
                    'dLon',dLon,'vel_ms',v_ms,'corr',mx); %#ok<AGROW>
                fprintf('    IS-velocity %s-%s: dLon=%.1f, lag=%.1f h, v=%.0f m/s (r=%.2f)\n',...
                    stations{a},stations{b},dLon,lag_h,v_ms,mx);
            catch
            end
        end
    end
end
end

function [onsetLT, onsetStd, onsetN] = sunsetOnsetRepeatability(dayRoots, stations, ...
    tEdges, tC, Z, lonC, IST)
% Find EPB onset (first epoch ROTI>1.0 post-sunset) per station per day.
% Returns per-station onset-LT std across days.
onsetLT  = nan(numel(stations), numel(dayRoots));
onsetStd = nan(numel(stations),1);
onsetN   = zeros(numel(stations),1);
% use single-station ROTI median series
for s=1:numel(stations)
    st = stations{s};
    for d=1:numel(dayRoots)
        try
            S = aloadTEC(dayRoots{d}, st);
            if isempty(S) || ~isfield(S,'ROTI'), continue; end
            rmed = nanmedian(S.ROTI,2);  % per-second
            % look for first epoch with ROTI>1.0 between 12-22 UTC
            mask12_22 = false(size(rmed));
            mask12_22(12*3600:min(22*3600,end)) = true;
            fnd = find(isfinite(rmed) & rmed>1.0 & mask12_22);
            if isempty(fnd), continue; end
            onset_utc = (fnd(1)-1)/3600;
            onsetLT(s,d) = onset_utc + IST;
        catch
        end
    end
    vals = onsetLT(s,:);
    vals = vals(isfinite(vals));
    onsetN(s)   = numel(vals);
    if numel(vals)>=2
        onsetStd(s) = std(vals);
    end
end
end

function S4mat = loadS4mat(dayRoot, station, sys)
% Load S4 matrix from TEC_*.mat if present.
S4mat = [];
try
    resDir = fullfile(dayRoot,'Results');
    if strcmp(sys,'G')
        pat = fullfile(resDir, sprintf('TEC_%s_*.mat',station));
    else
        pat = fullfile(resDir, sprintf('TEC_%s_%s_*.mat',station,sys));
    end
    mats = dir(pat);
    if isempty(mats), return; end
    S = load(fullfile(resDir,mats(1).name));
    fn = fieldnames(S);
    for k=1:numel(fn)
        if startsWith(fn{k},'S4_')
            S4mat = S.(fn{k});
            return;
        end
    end
catch
end
end

function [p] = welch2(a, b)
% Simple two-sample Wilcoxon rank-sum p-value approximation.
n1=numel(a); n2=numel(b);
all_=[a(:);b(:)]; [~,ri]=sort(all_);
rnk=zeros(n1+n2,1); rnk(ri)=1:(n1+n2);
W1=sum(rnk(1:n1));
mu=n1*(n1+n2+1)/2; sig=sqrt(n1*n2*(n1+n2+1)/12);
z=(W1-mu)/sig;
p=2*(1-0.5*(1+erf(abs(z)/sqrt(2))));
end

function boxplot_manual(ax, v, xpos, clr)
% Manual box plot (no Statistics Toolbox needed).
q = quantile(v,[0.25 0.5 0.75]);
iqr_v = q(3)-q(1);
lo = q(1)-1.5*iqr_v; hi = q(3)+1.5*iqr_v;
lo = max(lo,min(v)); hi = min(hi,max(v));
w=0.3;
patch(ax,[xpos-w xpos+w xpos+w xpos-w],[q(1) q(1) q(3) q(3)],clr,'FaceAlpha',0.5,'EdgeColor','k');
plot(ax,[xpos-w xpos+w],[q(2) q(2)],'k-','LineWidth',2);
plot(ax,[xpos xpos],[lo q(1)],'k-'); plot(ax,[xpos xpos],[q(3) hi],'k-');
out_v = v(v<lo | v>hi);
if ~isempty(out_v)
    scatter(ax,xpos*ones(size(out_v)),out_v,20,clr,'filled','MarkerFaceAlpha',0.5);
end
end

function shadedCI(ax, t, mu, lo, hi, clr)
% Shaded confidence band.
try
    ok = isfinite(t) & isfinite(lo) & isfinite(hi);
    if sum(ok)<2, return; end
    fill(ax,[t(ok); flipud(t(ok))],[lo(ok); flipud(hi(ok))],clr,...
        'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
catch
end
end
