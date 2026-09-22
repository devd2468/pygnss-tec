function aplotIPPmap(dayRoot, station4, sysChar, utcHours, field, outDir)
% APLOTIPPMAP  True 2D interpolated TEC/ROTI maps from per-IPP geometry.
% Replaces station-level approximations with thin-shell IPP lat/lon from
% computeIPP.m (350 km shell) + natural-neighbor gridding.
%
% Inputs:
%   dayRoot   day folder, station4 4-char code, sysChar 'G'/'E'/'C'/'J'
%   utcHours  vector of snapshot hours, e.g. [14 16 18]
%   field     'vtec' (default) or 'roti'
%   outDir    default <dayRoot>/Results
    if nargin < 5 || isempty(field), field = 'vtec'; end
    if nargin < 6 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    ippFile = fullfile(outDir, sprintf('IPP_%s_%s.mat', station4, sysChar));
    if ~exist(ippFile, 'file')
        P = computeIPP(dayRoot, station4, sysChar, outDir);
    else
        P = load(ippFile);
    end
    if isempty(P) || ~isfield(P, 'lat')
        fprintf('aplotIPPmap: no IPP data for %s/%s\n', station4, sysChar);
        return;
    end
    if ~isfield(P, field)
        fprintf('aplotIPPmap: field %s missing, using vtec\n', field);
        field = 'vtec';
    end
    Z = P.(field);
    nP = numel(utcHours);
    fig = figure('Name', 'IPP map', ...
        'Position', [100 50 420*nP + 120, 640], 'Visible', 'off');
    for k = 1:nP
        hh = utcHours(k);
        rows = max(1, round(hh*3600)+(-60:60));  % +-60 s window
        rows = rows(rows >= 1 & rows <= size(Z, 1));
        xs = []; ys = []; vs = [];
        for c = 1:size(Z, 2)
            vv = Z(rows, c);
            ok = isfinite(vv) & isfinite(P.lat(rows, c)) & isfinite(P.lon(rows, c));
            if any(ok)
                xs = [xs; P.lon(rows(ok), c)]; %#ok<AGROW>
                ys = [ys; P.lat(rows(ok), c)]; %#ok<AGROW>
                vs = [vs; vv(ok)]; %#ok<AGROW>
            end
        end
        ax = subplot(1, nP, k); hold on; grid on;
        if numel(vs) >= 10
            [lonG, latG] = meshgrid(linspace(min(xs)-1, max(xs)+1, 80), ...
                                    linspace(min(ys)-1, max(ys)+1, 80));
            F = griddata(xs, ys, vs, lonG, latG, 'natural');
            pcolor(lonG, latG, F);
            shading interp;
        end
        if ~isempty(vs)
            scatter(xs, ys, 18, vs, 'filled', 'MarkerEdgeColor', 'k');
        end
        plot(P.reflon, P.reflat, 'k^', 'MarkerSize', 12, 'LineWidth', 2);
        text(P.reflon+0.3, P.reflat, station4, 'FontSize', 9, 'FontWeight', 'bold');
        % magnetic equator (dip ~0): EIA crest guides at +/-17.5 deg approx
        xl = xlim;
        plot(xl, [0 0], 'k-', 'LineWidth', 1.2);
        xlabel('Longitude'); ylabel('Latitude');
        title(sprintf('%02d:00 UTC  %s %s (n=%d IPP)', round(hh), station4, sysChar, numel(vs)));
        colorbar;
        if strcmp(field, 'vtec'), caxis([0 80]); else, caxis([0 1]); end
    end
    sgtitle(sprintf('%s %s true-IPP %s maps (350 km shell)', station4, sysChar, upper(field)));
    out = fullfile(outDir, sprintf('aplotIPPmap_%s_%s_%s.png', station4, sysChar, field));
    saveas(gcf, out);
    close;
    fprintf('aplotIPPmap saved: %s\n', out);
end
