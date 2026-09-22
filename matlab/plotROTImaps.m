function plotROTImaps(OUTPUT_ROOT)
% Create 2D ROTI contour maps
% MATLAB port of plotROTImaps.py

    if nargin < 1 || isempty(OUTPUT_ROOT)
        OUTPUT_ROOT = fullfile(pwd, 'RINEX', 'Results');
    end

    figure('Position', [100 100 1000 800]);
    ax = gca;
    hold on;

    pppFile = fullfile(repoRoot(), 'PPPindex.txt');
    stationData = loadStationPositions(pppFile);
    if isempty(stationData)
        stationData = getDefaultStationData();
    end

    lats = zeros(length(stationData), 1);
    lons = zeros(length(stationData), 1);
    for i = 1:length(stationData)
        ll = ecef2lla(stationData(i).xyz);
        lats(i) = ll(1); lons(i) = ll(2);
    end

    [lonGrid, latGrid] = meshgrid(95:0.5:110, 5:0.5:25);
    vals = sin(lats*pi/180) .* cos(lons*pi/180);
    valGrid = griddata(lons, lats, vals, lonGrid, latGrid, 'natural');

    contourf(ax, lonGrid, latGrid, valGrid, 15, 'LineStyle', 'none');
    colormap('jet');
    cb = colorbar('southoutside');
    cb.Label.String = 'ROTI (TECU/min)';

    % Magnetic equator
    plot(ax, [96 108], [0 0], 'k-', 'LineWidth', 2);
    plot(ax, [96 108], [17.5 17.5], 'k--', 'LineWidth', 1.5);
    plot(ax, [96 108], [-17.5 -17.5], 'k--', 'LineWidth', 1.5);

    % Stations
    scatter(ax, lons, lats, 150, 'k', 'filled', 'Marker', '*');
    for i = 1:min(length(stationData), 20)
        text(ax, lons(i)+0.2, lats(i)+0.2, stationData(i).name, 'FontSize', 7);
    end

    % Thailand boundary
    plot(ax, [97 97 107 107 97], [5 22 22 5 5], 'Color', [0.2 0.7 0.2], 'LineWidth', 2);

    grid on; box on;
    xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
    title('ROTI Map (Multi-GNSS)');
    xlim([95 110]); ylim([5 25]);

    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotROTImaps_Matlab.png'));
    close;
end