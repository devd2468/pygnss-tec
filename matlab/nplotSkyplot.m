function nplotSkyplot(dayRoot, station4, outDir, sysList)
% NPLOTSKYPLOT  Polar skyplots of IPP tracks colored by ROTI (one panel
% per constellation). theta = azimuth, rho = 90 - elevation.
    if nargin < 4 || isempty(sysList), sysList = {'G', 'E'}; end
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fig = figure('Name', 'skyplot', 'Position', [100 50 420*numel(sysList)+80, 520], 'Visible', 'off');
    for t = 1:numel(sysList)
        ax = subplot(1, numel(sysList), t, polaraxes);
        hold on;
        try, P = getIPP(dayRoot, station4, sysList{t});
        catch, P = [];
        end
        if ~isempty(P) && isfield(P, 'azim')
            th = deg2rad(P.azim(:));
            rh = 90 - P.elev(:);
            vv = P.roti(:);
            ok = isfinite(th) & isfinite(rh) & isfinite(vv) & rh <= 60;
            if any(ok)
                scatter(ax, th(ok), rh(ok), 6, vv(ok), 'filled', 'MarkerFaceAlpha', 0.6);
            end
        end
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        ax.RLim = [0 60];
        ax.RTick = [0 15 30 45 60];
        ax.RTickLabel = {'90', '75', '60', '45', '30'};
        colormap(ax, 'jet');
        caxis([0 1]);
        title(sprintf('%s %s (elev>30\\circ)', station4, sysList{t}));
    end
    colorbar('eastoutside');
    sgtitle(sprintf('%s IPP sky tracks colored by ROTI', station4));
    out = fullfile(outDir, sprintf('nplotSkyplot_%s.png', station4));
    figsavesafe(gcf, out);
    fprintf('nplotSkyplot saved: %s\n', out);
end
