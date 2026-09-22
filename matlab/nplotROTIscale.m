function nplotROTIscale(dayRoots, stations, outDir)
% NPLOTROTISCALE  Vankadara-style 4-level ROTI occurrence climatology:
% quiet [0,.25), weak [.25,.5), moderate [.5,1), strong [1,inf),
% percentage of valid ROTI samples per station-day (+ pooled row).
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    edges = [0, 0.25, 0.5, 1.0, Inf];
    P = nan(numel(stations), numel(dayRoots), 4);
    for d = 1:numel(dayRoots)
        for s = 1:numel(stations)
            S = aloadTEC(dayRoots{d}, stations{s});
            if isempty(S) || ~isfield(S, 'ROTI'), continue; end
            v = S.ROTI(:);
            v = v(isfinite(v));
            if isempty(v), continue; end
            for k = 1:4
                P(s, d, k) = 100 * sum(v >= edges(k) & v < edges(k+1)) / numel(v);
            end
        end
    end
    Pm = nanmean(P, 2);  % pooled over days -> stations x 4
    fig = figure('Name', 'rotiscale', 'Position', [100 50 1100 600], 'Visible', 'off');
    subplot(1, 2, 1);
    bar(squeeze(nanmean(P, 2))', 'stacked');
    set(gca, 'XTickLabel', stations, 'XTickLabelRotation', 30);
    ylabel('% of ROTI samples');
    title('ROTI occurrence by station (5-day pooled)');
    legend({'quiet <0.25', 'weak 0.25-0.5', 'moderate 0.5-1', 'strong >=1'}, 'Location', 'best');
    grid on;
    subplot(1, 2, 2);
    imagesc(1:numel(dayRoots), 1:numel(stations), squeeze(nanmean(P(:, :, 3:4), 3)));
    set(gca, 'YTick', 1:numel(stations), 'YTickLabel', stations);
    set(gca, 'XTick', 1:numel(dayRoots), 'XTickLabel', ...
        arrayfun(@(d) sprintf('D%d', d), 1:numel(dayRoots), 'UniformOutput', false));
    xlabel('Day'); ylabel('Station');
    title('% moderate+strong ROTI (>=0.5)');
    colorbar;
    out = fullfile(outDir, 'nplotROTIscale.png');
    figsavesafe(gcf, out);
    % CSV table for the paper
    fid = fopen(fullfile(outDir, 'roti_scale_table.csv'), 'w');
    fprintf(fid, 'station,quiet,weak,moderate,strong\n');
    for s = 1:numel(stations)
        fprintf(fid, '%s,%.2f,%.2f,%.2f,%.2f\n', stations{s}, Pm(s, 1), Pm(s, 2), Pm(s, 3), Pm(s, 4));
    end
    fclose(fid);
    fprintf('nplotROTIscale saved + roti_scale_table.csv\n');
end
