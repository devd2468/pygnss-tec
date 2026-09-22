function aplotS4(dayRoots, stations, outDir, varargin)
%{
 ================================================================================
  aplotS4  —  3-panel S4 amplitude-scintillation article figure
 ================================================================================

 Produces a publication-quality 3-panel figure comparing S4 (amplitude
 scintillation) with ROTI (phase scintillation) for one or more station-days.

 Panel 1 — S4 + ROTI dual-axis time series
   Left  axis : S4 index (60-s window, dimensionless) — amplitude scintillation
   Right axis : ROTI median (TECU/min)              — phase scintillation
   Horizontal dashed lines: s4Weak=0.15, s4Moderate=0.30, s4Strong=0.50

 Panel 2 — S4 vs ROTI scatter coloured by elevation angle
   Reveals the amplitude-phase relationship and elevation dependence.
   Each point = one (PRN, epoch) pair where both are finite.
   Colour: elevation [deg], colormap 'jet' (low elev = blue, high = red).

 Panel 3 — S4 occurrence rate by station vs local solar time
   Histogram of S4 > s4Weak epochs binned by local solar time (LST = UTC +
   station_lon_hours). Bars coloured per station. Highlights post-sunset
   enhancement typical of EPB (equatorial plasma bubble) events.

 Inputs
 ------
   dayRoots  - cell array of day-folder paths, each containing Results/
               TEC_<ST>_*.mat with S4_<ST>_* variable.
               Example: {'D:\IGS\2025-05-02', 'D:\IGS\2025-05-04'}
   stations  - cell array of 4-char station codes, e.g. {'BHPL','LCK4'}
               Use {} to load all stations found in each dayRoot/Results/.
   outDir    - output folder for saved PNG (created if absent).
   varargin  - optional name/value pairs:
               'titleStr'  char label for sgtitle (default: date from file)
               'figWidth'  figure width  px (default 1400)
               'figHeight' figure height px (default 420 per panel = 1260)
               'saveFig'   true|false (default true)

 Output
 ------
   aplotS4_<YYYYMMDD>.png  saved to outDir
   Figure handle returned (if nargout > 0)

 Example
 -------
   aplotS4({'D:\IGS\2025-05-04'}, {'BHPL'}, 'D:\IGS\ArticleFigs')
   aplotS4({'D:\IGS\2025-05-02','D:\IGS\2025-05-04'}, {}, outDir, ...
           'titleStr','May 2025 EPB Night','saveFig',true)

 Dependencies
 ------------
   figsavesafe.m  (safe figure export wrapper)
   config.m       (reads s4Weak/s4Moderate/s4Strong thresholds)
   PPPindex.txt   (station longitudes for LST calculation)

 CSSRG Laboratory, KMITL | aplotS4 — B1 (S4) block, Sep 2026
%}

%% ---------- defaults & option parsing ----------
p_titleStr  = '';
p_figWidth  = 1400;
p_figHeight = 1260;
p_saveFig   = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'titlestr',  p_titleStr  = varargin{k+1};
        case 'figwidth',  p_figWidth  = varargin{k+1};
        case 'figheight', p_figHeight = varargin{k+1};
        case 'savefig',   p_saveFig   = varargin{k+1};
    end
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
cfg = config();   % read thresholds from single config source

%% ---------- load S4 / ROTI / prm across dayRoots x stations ----------
S4_all   = [];   % [N x 1] S4 values
ROTI_all = [];   % [N x 1] ROTI values (matched)
ELEV_all = [];   % [N x 1] elevation [deg]
LST_all  = [];   % [N x 1] local solar time [h]
STAT_all = {};   % {N x 1} station labels

% Station longitude lookup (for LST)
try
    sp = loadStationPositions(fullfile(repoRoot(), 'PPPindex.txt'));
catch
    sp = getDefaultStationData();
end

dateStr = '';  % for auto-title

for di = 1:numel(dayRoots)
    resDir = fullfile(dayRoots{di}, 'Results');
    if ~exist(resDir, 'dir')
        fprintf('aplotS4: Results dir not found: %s\n', resDir);
        continue;
    end

    % Determine stations to load
    if isempty(stations)
        mats = dir(fullfile(resDir, 'TEC_*_????_??_??.mat'));
        stList = {};
        for k = 1:numel(mats)
            tok = regexp(mats(k).name, '^TEC_([A-Z0-9]{4})_\d{4}', 'tokens', 'once');
            if ~isempty(tok), stList{end+1} = tok{1}; end %#ok<AGROW>
        end
        stList = unique(stList);
    else
        stList = stations;
    end

    for si = 1:numel(stList)
        st = stList{si};
        % Find the .mat file with S4 (saved by TECcalculationRINEX304_OEM7)
        pat = fullfile(resDir, sprintf('TEC_%s_*.mat', st));
        mats = dir(pat);
        if isempty(mats)
            fprintf('aplotS4: No TEC mat for %s in %s\n', st, resDir);
            continue;
        end
        fname = fullfile(resDir, mats(1).name);

        % Date from filename for title
        tok = regexp(mats(1).name, '_(\d{4})_(\d{2})_(\d{2})\.mat$', 'tokens', 'once');
        if ~isempty(tok) && isempty(dateStr)
            dateStr = sprintf('%s-%s-%s', tok{1}, tok{2}, tok{3});
        end

        % Load variables from .mat
        try
            S = load(fname);
        catch ME
            fprintf('aplotS4: load failed %s: %s\n', fname, ME.message);
            continue;
        end

        % --- find S4, ROTI, prm variable names (year_month_date suffix) ---
        fnames_mat = fieldnames(S);
        s4_name   = findVar(fnames_mat, sprintf('S4_%s',   st));
        roti_name  = findVar(fnames_mat, sprintf('ROTI_%s', st));
        prm_name   = findVar(fnames_mat, sprintf('prm_%s',  st));

        % Fallback: no-station-prefix (GPS engine saves without station tag)
        if isempty(s4_name)
            s4_name = findVarAny(fnames_mat, 'S4_');
        end
        if isempty(roti_name)
            roti_name = findVarAny(fnames_mat, 'ROTI_');
        end
        if isempty(prm_name)
            prm_name = findVarAny(fnames_mat, 'prm_');
        end

        if isempty(s4_name)
            fprintf('aplotS4: S4 variable not found in %s (run pipeline first)\n', fname);
            continue;
        end

        S4_mat   = S.(s4_name);                           % 86400 x NPRN
        ROTI_mat = [];
        prm_mat  = struct('elevation', nan(86400, size(S4_mat, 2)));
        if ~isempty(roti_name),  ROTI_mat = S.(roti_name); end
        if ~isempty(prm_name),   prm_mat  = S.(prm_name);  end

        % Station longitude for LST
        lon_deg = 0;
        ix = find(strcmp({sp.name}, st), 1);
        if isempty(ix), ix = find(strncmp({sp.name}, st, 4), 1); end
        if ~isempty(ix)
            lla = ecef2lla(sp(ix).xyz);
            lon_deg = lla(2);
        end

        % Flatten matrices: collect valid (S4, ROTI, elev, LST) tuples
        [nT, nPRN] = size(S4_mat);
        for prn = 1:nPRN
            s4col   = S4_mat(:, prn);
            rotcol  = nan(nT, 1);
            elevcol = nan(nT, 1);
            if ~isempty(ROTI_mat) && size(ROTI_mat, 2) >= prn
                rotcol = ROTI_mat(:, prn);
            end
            if isfield(prm_mat, 'elevation') && size(prm_mat.elevation, 2) >= prn
                elevcol = prm_mat.elevation(:, prn);
            end

            % Filter: both S4 and elevation must be finite
            ok = isfinite(s4col) & isfinite(elevcol) & (elevcol > 0);
            if ~any(ok), continue; end

            t_sec = (find(ok) - 1);   % 0-based seconds from midnight
            lst   = mod(t_sec / 3600 + lon_deg / 15, 24);  % local solar time [h]

            S4_all   = [S4_all;   s4col(ok)];   %#ok<AGROW>
            if ~isempty(ROTI_mat)
                ROTI_all = [ROTI_all; rotcol(ok)]; %#ok<AGROW>
            else
                ROTI_all = [ROTI_all; nan(sum(ok), 1)]; %#ok<AGROW>
            end
            ELEV_all = [ELEV_all;  elevcol(ok)]; %#ok<AGROW>
            LST_all  = [LST_all;   lst(:)];       %#ok<AGROW>
            STAT_all = [STAT_all;  repmat({st}, sum(ok), 1)]; %#ok<AGROW>
        end
    end
end

if isempty(S4_all)
    warning('aplotS4: no valid S4 data found; nothing to plot.');
    return;
end

%% ---------- build figure ----------
fig = figure('Name', 'S4 Amplitude Scintillation', ...
    'Position', [100 50 p_figWidth p_figHeight], 'Visible', 'off');

stUniq = unique(STAT_all);
nSt    = numel(stUniq);
cmap_st = lines(max(nSt, 1));

%% --- PANEL 1: S4 + ROTI time series (dual Y-axis, per station) ----------
ax1 = subplot(3, 1, 1);
yyaxis left
hold on; grid on;
for si = 1:nSt
    st = stUniq{si};
    mask = strcmp(STAT_all, st);
    % Use median across PRNs per second (aggregate to time series)
    t_idx = LST_all(mask);
    s4v   = S4_all(mask);
    % Bin into 1-min UTC buckets, take median
    tbins = 0:1/60:24;
    s4_med = nan(numel(tbins)-1, 1);
    for b = 1:numel(tbins)-1
        inb = t_idx >= tbins(b) & t_idx < tbins(b+1);
        if any(inb), s4_med(b) = nanmedian(s4v(inb)); end
    end
    t_center = (tbins(1:end-1) + tbins(2:end)) / 2;
    plot(t_center, s4_med, 'LineWidth', 1.2, 'Color', cmap_st(si, :), ...
        'DisplayName', st);
end
yline(cfg.s4Weak,     'b--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(cfg.s4Moderate, 'r--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(cfg.s4Strong,   'k-',  'LineWidth', 1.0, 'HandleVisibility', 'off');
ylabel('S4 (dim-less)');
ylim([0, max(1.0, max(S4_all)*1.1)]);
set(ax1, 'YColor', 'k');

yyaxis right
if any(isfinite(ROTI_all))
    for si = 1:nSt
        st = stUniq{si};
        mask = strcmp(STAT_all, st);
        t_idx = LST_all(mask);
        rv    = ROTI_all(mask);
        tbins = 0:1/60:24;
        r_med = nan(numel(tbins)-1, 1);
        for b = 1:numel(tbins)-1
            inb = t_idx >= tbins(b) & t_idx < tbins(b+1);
            if any(inb), r_med(b) = nanmedian(rv(inb)); end
        end
        t_center = (tbins(1:end-1) + tbins(2:end)) / 2;
        plot(t_center, r_med, '--', 'LineWidth', 0.9, 'Color', cmap_st(si, :) * 0.7, ...
            'HandleVisibility', 'off');
    end
    ylabel('ROTI (TECU/min)');
    yline(cfg.rotiMinor,    'b:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
    yline(cfg.rotiModerate, 'r:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
    yline(cfg.rotiSevere,   'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
xlim([0 24]);
xlabel('Local Solar Time (h)');
title('S4 (solid) and ROTI (dashed) time series — amplitude vs phase scintillation');
legend(stUniq, 'Location', 'northwest');
set(ax1, 'XTick', 0:2:24);

%% --- PANEL 2: S4 vs ROTI scatter coloured by elevation -----------------
subplot(3, 1, 2);
hold on; grid on;
if any(isfinite(ROTI_all))
    ok2 = isfinite(S4_all) & isfinite(ROTI_all) & isfinite(ELEV_all);
    if any(ok2)
        sc = scatter(ROTI_all(ok2), S4_all(ok2), 4, ELEV_all(ok2), ...
            'filled', 'MarkerFaceAlpha', 0.4);
        colormap(gca, jet(256));
        cb = colorbar;
        cb.Label.String = 'Elevation (deg)';
        caxis([0, 90]);
        % Threshold lines
        xline(cfg.rotiMinor,    'b--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        xline(cfg.rotiModerate, 'r--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        xline(cfg.rotiSevere,   'k-',  'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(cfg.s4Weak,       'b--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(cfg.s4Moderate,   'r--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(cfg.s4Strong,     'k-',  'LineWidth', 0.8, 'HandleVisibility', 'off');
        % Correlation annotation
        rho = corr(ROTI_all(ok2), S4_all(ok2), 'rows', 'complete');
        text(0.98, 0.97, sprintf('r = %.3f', rho), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 9);
    end
else
    text(0.5, 0.5, 'ROTI not available', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
end
xlabel('ROTI (TECU/min)');
ylabel('S4 (dim-less)');
title('S4 vs ROTI scatter (coloured by elevation) — the key science plot');

%% --- PANEL 3: S4 occurrence by station vs local solar time -------------
subplot(3, 1, 3);
hold on; grid on;
edges  = 0:0.5:24;   % 30-min LST bins
for si = 1:nSt
    st = stUniq{si};
    mask = strcmp(STAT_all, st) & (S4_all > cfg.s4Weak);
    if ~any(mask), continue; end
    lst_ev = LST_all(mask);
    N = histcounts(lst_ev, edges);
    % Normalise to fractional occurrence rate per bin
    N_tot = histcounts(LST_all(strcmp(STAT_all, st)), edges);
    N_tot(N_tot == 0) = 1;          % avoid /0
    rate = N ./ N_tot;
    t_c = (edges(1:end-1) + edges(2:end)) / 2;
    bar(t_c, rate, 0.9, 'FaceColor', cmap_st(si, :), ...
        'FaceAlpha', 0.65, 'DisplayName', st, 'EdgeColor', 'none');
end
xlim([0 24]);
xlabel('Local Solar Time (h)');
ylabel('S4 > 0.15 occurrence rate');
title(sprintf('S4 occurrence by station vs local time (threshold %.2f)', cfg.s4Weak));
legend(stUniq, 'Location', 'northwest');
set(gca, 'XTick', 0:2:24);

%% ---------- super-title & save ----------
if isempty(p_titleStr)
    if ~isempty(dateStr)
        p_titleStr = sprintf('S4 Amplitude Scintillation — %s', dateStr);
    else
        p_titleStr = 'S4 Amplitude Scintillation';
    end
end
sgtitle(p_titleStr, 'FontSize', 13, 'FontWeight', 'bold');

if p_saveFig
    if ~isempty(dateStr)
        outName = fullfile(outDir, sprintf('aplotS4_%s.png', strrep(dateStr,'-','')));
    else
        outName = fullfile(outDir, 'aplotS4.png');
    end
    figsavesafe(fig, outName);
    fprintf('aplotS4: saved %s\n', outName);
end

if nargout < 1
    close(fig);
end

end  % aplotS4

%% ============================== helpers =================================
function vname = findVar(fnames, prefix)
% Return first fieldname that starts with 'prefix_'
vname = '';
for k = 1:numel(fnames)
    if strncmp(fnames{k}, [prefix '_'], length(prefix)+1)
        vname = fnames{k};
        return;
    end
end
end

function vname = findVarAny(fnames, prefix)
% Return first fieldname that starts with 'prefix' (no underscore constraint)
vname = '';
for k = 1:numel(fnames)
    if strncmp(fnames{k}, prefix, length(prefix))
        vname = fnames{k};
        return;
    end
end
end
