function plotMappingRange(ax, stationLon, stationLat, radius)
% Plot range circle around station
% MATLAB port of plotMappingRange.py

    if nargin < 4 || isempty(radius)
        radius = 500; % km
    end

    R_earth = 6371; % km
    delta = radius / R_earth; % degrees angular distance

    angles = linspace(0, 2*pi, 100);
    circleLat = stationLat + delta * sind(angles);
    circleLon = stationLon + delta * sind(angles) ./ cosd(stationLat);

    if isempty(ax)
        figure('Position', [100 100 800 600]);
        ax = gca;
    end

    plot(ax, circleLon, circleLat, 'k-', 'LineWidth', 1.5);
    hold on;
    plot(ax, stationLon, stationLat, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    text(ax, stationLon + 0.3, stationLat, ...
        sprintf('Station (%d km radius)', radius), 'FontSize', 10);
    grid on;
    xlabel('Longitude'); ylabel('Latitude');
    title('Station Observation Range');
end