function aplotIPPgain(dayRoot, station4, outDir)
% APLOTIPPGAIN  IPP-count gain from adding Galileo/BeiDou/QZSS (Phase 4
% headline figure). Counts satellites with elevation > 30 deg per epoch
% for GPS alone vs GPS+E vs GPS+E+C vs GPS+E+C+J.
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    base = dir(fullfile(dayRoot, 'Results', sprintf('MultiGNSS_%s_2*.mat', station4)));
    if isempty(base)
        base = dir(fullfile(dayRoot, 'Results', sprintf('MultiGNSS_%s_*.mat', station4)));
    end
    if isempty(base)
        fprintf('aplotIPPgain: no data for %s in %s\n', station4, dayRoot);
        return;
    end
    % GPS base file = shortest name match without extra _SYS_ tag
    gFile = '';
    for k = 1:numel(base)
        if isempty(regexp(base(k).name, 'MultiGNSS_[A-Z0-9]+_[ERCJ]_', 'once'))
            gFile = fullfile(base(k).folder, base(k).name);
            break;
        end
    end
    if isempty(gFile)
        fprintf('aplotIPPgain: no GPS base file for %s\n', station4);
        return;
    end
    G = load(gFile);
    if ~isfield(G, 'prm') || ~isfield(G.prm, 'elevation')
        fprintf('aplotIPPgain: no elevation in %s\n', gFile);
        return;
    end
    n = size(G.prm.elevation, 1);
    cnt = zeros(n, 4);  % G, +E, +E+C, +E+C+J
    cnt(:, 1) = sum(G.prm.elevation > 30, 2);
    tags = {'E', 'C', 'J'};
    acc = cnt(:, 1);
    for t = 1:numel(tags)
        pat = sprintf('MultiGNSS_%s_%s_*.mat', station4, tags{t});
        d = dir(fullfile(dayRoot, 'Results', pat));
        add = zeros(n, 1);
        for k = 1:numel(d)
            try
                S = load(fullfile(d(k).folder, d(k).name));
                if isfield(S, 'prm') && isfield(S.prm, 'elevation')
                    E = S.prm.elevation;
                    if size(E, 1) ~= n
                        % resample to GPS grid length
                        Eq = nan(n, size(E, 2));
                        m = min(n, size(E, 1));
                        Eq(1:m, :) = E(1:m, :);
                        E = Eq;
                    end
                    add = add + sum(E > 30, 2);
                end
            catch
            end
        end
        acc = acc + add;
        cnt(:, t+1) = acc;
    end
    tH = (0:n-1) / 3600;
    fig = figure('Name', 'IPP gain', 'Position', [100 50 1100 650], 'Visible', 'off');
    hold on; grid on;
    area(tH, cnt(:, 4), 'FaceColor', [0.75 0.85 1], 'EdgeColor', 'none');
    area(tH, cnt(:, 3), 'FaceColor', [0.65 0.9 0.65], 'EdgeColor', 'none');
    area(tH, cnt(:, 2), 'FaceColor', [1 0.85 0.5], 'EdgeColor', 'none');
    plot(tH, cnt(:, 1), 'b-', 'LineWidth', 1.6);
    xlabel('Universal Time (hour)');
    ylabel('Satellites elev > 30° (IPP count)');
    title(sprintf('%s  IPP density: GPS vs +Galileo/BeiDou/QZSS', station4));
    legend('GPS+E+C+J', 'GPS+E+C', 'GPS+E', 'GPS only', 'Location', 'best');
    xlim([tH(1), tH(end)]);
    g0 = nanmean(cnt(:, 1)); g3 = nanmean(cnt(:, 4));
    text(0.02, 0.92, sprintf('mean GPS-only: %.1f  |  mean +ECJ: %.1f  (x%.2f)', ...
        g0, g3, g3/max(g0, eps)), 'Units', 'normalized', 'FontSize', 11);
    out = fullfile(outDir, sprintf('aplotIPPgain_%s.png', station4));
    saveas(gcf, out);
    close;
    fprintf('aplotIPPgain saved: %s (gain x%.2f)\n', out, g3/max(g0, eps));
end
