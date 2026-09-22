function plotScatterAndInterpolation(OUTPUT_ROOT)
% IPP scatter with spatial interpolation
% MATLAB port of plotScatterAndInterpolation.py

    if nargin < 1 || isempty(OUTPUT_ROOT)
        OUTPUT_ROOT = fullfile(pwd, 'RINEX', 'Results');
    end

    figure('Position', [100 100 1400 600]);

    matFiles = dir(fullfile(OUTPUT_ROOT, 'MultiGNSS_*.mat'));
    if isempty(matFiles)
        matFiles = dir(fullfile(OUTPUT_ROOT, 'TEC_*.mat'));
    end
    if isempty(matFiles)
        fprintf('No data files found.\n');
        return;
    end

    data = load(fullfile(matFiles(1).folder, matFiles(1).name));
    if ~isfield(data, 'TEC')
        fprintf('No TEC data.\n');
        return;
    end

    % Left: STEC scatter
    subplot(1, 2, 1);
    if isfield(data.TEC, 'STECp')
        stecl_data = data.TEC.STECp;
    else
        stecl_data = data.TEC.STECl;
    end
    if isfield(data, 'prm') && isfield(data.prm, 'elevation')
        elev_data = data.prm.elevation;
    else
        elev_data = nan(size(stecl_data));
    end

    for prn = 1:32
        valid = ~isnan(stecl_data(:, prn)) & ~isnan(elev_data(:, prn)) & (elev_data(:, prn) > 30);
        if any(valid)
            scatter(stecl_data(valid, prn), elev_data(valid, prn), 15, ...
                stecl_data(valid, prn), 'filled', 'MarkerEdgeColor', 'none');
            hold on;
        end
    end
    colorbar;
    xlabel('STEC (TECU)'); ylabel('Elevation (°)');
    title('STEC vs Elevation (Elev > 30°)'); grid on;

    % Right: Interpolated grid
    subplot(1, 2, 2);
    stationData = loadStationPositions(fullfile(repoRoot(), 'PPPindex.txt'));
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

    contourf(lonGrid, latGrid, valGrid, 15, 'LineStyle', 'none');
    colormap('jet');
    colorbar; hold on;

    scatter(lons, lats, 100, 'k', 'filled', 'Marker', 's');
    for i = 1:min(length(stationData), 15)
        text(lons(i)+0.2, lats(i)+0.2, stationData(i).name, 'FontSize', 7);
    end

    grid on; xlabel('Longitude'); ylabel('Latitude');
    title('VTEC Spatial Interpolation (Grid)');

    sgtitle('IPP Points and Interpolation (Multi-GNSS)', 'FontSize', 14);

    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotScatterAndInterpolation_Matlab.png'));
    close;
end