function plotROTIdtrend(OUTPUT_ROOT)
% ROTI trend comparison with different smoothing windows
% MATLAB port of plotROTIdtrend.py

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

    data = load(fullfile(matFiles(1).folder, matFiles(1).name));
    if ~isfield(data, 'ROTI')
        fprintf('No ROTI data.\n');
        return;
    end

    [time_roti, roti_median] = roti5min(data.ROTI);

    figure('Position', [100 100 1400 700]);
    hold on;

    plot(time_roti, roti_median, 'k-', 'LineWidth', 1.5, ...
        'DisplayName', 'Original (5-min)');

    for delta_min = [5, 15, 30, 60]
        windowSize = max(1, round(delta_min / 5));
        smoothed = smooth(roti_median, windowSize, 'moving');
        plot(time_roti, smoothed, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%d min window', delta_min));
    end

    yline(0.2, 'b--', 'Minor');
    yline(0.5, 'r--', 'Moderate');
    yline(1.0, 'k-', 'Severe');

    xlabel('Time (UT)'); ylabel('ROTI (TECU/min)');
    title('ROTI Trend Comparison with Different Smoothing Windows');
    legend('Location', 'northeast');
    grid on; xlim([0 24]); ylim([0 3]);

    saveas(gcf, fullfile(OUTPUT_ROOT, 'plotROTIdtrend_Matlab.png'));
    close;
end