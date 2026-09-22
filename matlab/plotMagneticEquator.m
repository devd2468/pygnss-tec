function plotMagneticEquator(ax, dateStr, latLim, lonLim)
% Plots magnetic equator / EIA crest lines
% MATLAB port of plotMagneticEquator.py

    if nargin < 2 || isempty(dateStr)
        dateStr = '2021';
    end
    if nargin < 3 || isempty(latLim)
        latLim = [-30, 30];
    end
    if nargin < 4 || isempty(lonLim)
        lonLim = [-100, 100];
    end

    if isempty(ax)
        figure('Position', [100 100 1000 600]);
        ax = gca;
    end

    lons = linspace(lonLim(1), lonLim(2), 360);

    % Equator
    plot(ax, lons, zeros(size(lons)), 'k-', 'LineWidth', 1.5);
    hold on;

    % EIA crests at ~±17.5° latitude (standard approximation for Thailand region)
    plot(ax, lons, 17.5*ones(size(lons)), '--', 'LineWidth', 2, ...
        'Color', [0.9 0.3 0.2], 'DisplayName', 'EIA Northern Crest');
    plot(ax, lons, -17.5*ones(size(lons)), '--', 'LineWidth', 2, ...
        'Color', [0.9 0.3 0.2], 'DisplayName', 'EIA Southern Crest');

    grid on;
    xlabel('Longitude (°E)');
    ylabel('Latitude (°N)');
    title(sprintf('Equatorial Ionization Anomaly Model (%s)', dateStr));
    legend('Location', 'southeast');
    ylim(latLim);

    % Fill EIA region
    fill([lons fliplr(lons)], [17.5*ones(size(lons)) fliplr(-17.5*ones(size(lons)))], ...
        [0.9 0.8 0.6], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    legend('Equator', 'EIA Northern Crest', 'EIA Southern Crest', 'EIA Region');
end