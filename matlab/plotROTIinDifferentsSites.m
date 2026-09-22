function plotROTIinDifferentsSites(OUTPUT_ROOT)
% Compare ROTI across all processed stations
% MATLAB port of plotROTIinDifferentsSites.py

    if nargin < 1 || isempty(OUTPUT_ROOT)
        OUTPUT_ROOT = fullfile(pwd, 'RINEX', 'Results');
    end

    matFiles = dir(fullfile(OUTPUT_ROOT, 'MultiGNSS_*.mat'));
    if isempty(matFiles)
        matFiles = dir(fullfile(OUTPUT_ROOT, 'TEC_*.mat'));
    end
    if isempty(matFiles)
        fprintf('No data files found.\n');
        return;
    end

    numStations = min(length(matFiles), 10);
    fig = figure('Position', [100 100 1000, 300 + 150 * numStations]);

    for s = 1:numStations
        data = load(fullfile(matFiles(s).folder, matFiles(s).name));

        tokens = regexp(matFiles(s).name, 'MultiGNSS_(\w+)_', 'tokens');
        if ~isempty(tokens)
            stationName = tokens{1}{1};
        else
            tokens = regexp(matFiles(s).name, 'TEC_(\w+)_', 'tokens');
            if ~isempty(tokens)
                stationName = tokens{1}{1};
            else
                stationName = sprintf('ST%02d', s);
            end
        end

        subplot(numStations, 1, s);

        if isfield(data, 'ROTI')
            [time_roti, roti_median] = roti5min(data.ROTI);

            plot(time_roti, roti_median, 'k.', 'MarkerSize', 3);
            hold on;
            yline(0.2, 'b--', 'Minor');
            yline(0.5, 'r--', 'Moderate');
            yline(1.0, 'k-', 'Severe');
        end

        ylabel('ROTI'); xlabel('Time (UTC)');
        title(sprintf('%s (Multi-GNSS)', stationName));
        xlim([0 24]); ylim([0 5]); grid on;
    end

    sgtitle('ROTI Time Series - All Stations', 'FontSize', 14);

    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotROTIinDifferentsSites_Matlab.png'));
    close;
end