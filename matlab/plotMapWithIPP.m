function plotMapWithIPP(OUTPUT_ROOT)
% Plot ROTI with IPP positions on map
% MATLAB port of plotMapWithIPP.py

    if nargin < 1 || isempty(OUTPUT_ROOT)
        OUTPUT_ROOT = fullfile(pwd, 'RINEX', 'Results');
    end

    figure('Position', [100 100 1000 800]);
    ax = gca;
    hold on;

    pppFile = fullfile(repoRoot(), 'PPPindex.txt');
    stationData = loadStationPositions(pppFile);

    if ~isempty(stationData)
        lats = zeros(length(stationData), 1);
        lons = zeros(length(stationData), 1);
        vals = zeros(length(stationData), 1);

        for i = 1:length(stationData)
            ll = ecef2lla(stationData(i).xyz);
            lats(i) = ll(1);
            lons(i) = ll(2);

            matFiles = dir(fullfile(OUTPUT_ROOT, 'MultiGNSS_*.mat'));
            if isempty(matFiles)
                matFiles = dir(fullfile(OUTPUT_ROOT, 'TEC_*.mat'));
            end
            if ~isempty(matFiles)
                try
                    data = load(fullfile(matFiles(1).folder, matFiles(1).name));
                    if isfield(data, 'TEC')
                        vals(i) = nanmean(data.TEC.vertical(:, i));
                    end
                catch
                    vals(i) = 0;
                end
            end
        end

        if any(vals ~= 0)
            scatter(ax, lons, lats, 300, vals, 'filled', ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8);
            cb = colorbar('southoutside');
            cb.Label.String = 'VTEC (TECU)';
            cb.Ticks = [0 20 40 60 80];
        else
            scatter(ax, lons, lats, 300, 'b', 'filled', 'MarkerEdgeColor', 'k');
        end

        for i = 1:length(stationData)
            text(ax, lons(i)+0.3, lats(i)+0.3, stationData(i).name, ...
                'FontSize', 8, 'FontWeight', 'bold');
        end
    end

    % Magnetic equator
    lons_range = linspace(96, 108, 100);
    plot(ax, lons_range, zeros(size(lons_range)), 'k-', 'LineWidth', 1.5);
    plot(ax, lons_range, 17.5*ones(size(lons_range)), 'r--', 'LineWidth', 1.5);
    plot(ax, lons_range, -17.5*ones(size(lons_range)), 'r--', 'LineWidth', 1.5);

    % Thailand boundary
    thailand_lat = [5, 22, 22, 5, 5];
    thailand_lon = [97, 97, 107, 107, 97];
    plot(ax, thailand_lon, thailand_lat, 'Color', [0.2 0.7 0.2], 'LineWidth', 2);

    grid on; box on;
    xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
    title('IPP Locations with VTEC Values (Multi-GNSS)');
    xlim([95, 110]); ylim([5, 25]);

    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotMapWithIPP_Matlab.png'));
    close;
end