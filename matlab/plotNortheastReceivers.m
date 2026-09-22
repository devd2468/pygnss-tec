function plotNortheastReceivers(OUTPUT_ROOT)
% Plot NE Thailand receiver network
% MATLAB port of plotNortheastReceivers.py

    if nargin < 1 || isempty(OUTPUT_ROOT)
        OUTPUT_ROOT = fullfile(pwd, 'RINEX', 'Results');
    end

    stationData = loadStationPositions(fullfile(repoRoot(), 'PPPindex.txt'));
    if isempty(stationData)
        stationData = getDefaultStationData();
    end

    ll = zeros(length(stationData), 3);
    for i = 1:length(stationData)
        ll(i, :) = ecef2lla(stationData(i).xyz);
    end

    lats = ll(:, 1);
    lons = ll(:, 2);
    names = {stationData.name};

    figure('Position', [100 100 1200 800]);
    ax = gca;
    hold on;

    % Stations as green circles with black edge
    scatter(ax, lons, lats, 200, [0.1 0.8 0.1], 'filled', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8);

    for i = 1:length(names)
        text(ax, lons(i)+0.2, lats(i), names{i}, ...
            'FontSize', 9, 'FontWeight', 'bold');
    end

    lat_min = min(lats) - 1; lat_max = max(lats) + 1;
    lon_min = min(lons) - 1; lon_max = max(lons) + 1;
    rectangle('Position', [lon_min lat_min lon_max-lon_min lat_max-lat_min], ...
        'EdgeColor', 'k', 'LineWidth', 4);

    grid on; box on;
    xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
    title('GNSS Receiver Network (NE Thailand) - Multi-GNSS', 'FontSize', 14);

    dummy = scatter(ax, [], [], 200, [0.1 0.8 0.1], 'filled', 'MarkerEdgeColor', 'k');
    legend(dummy, 'GNSS Receivers', 'Location', 'northwest');
    legend('boxoff');

    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotNortheastReceivers_Matlab.png'));
    close;
end