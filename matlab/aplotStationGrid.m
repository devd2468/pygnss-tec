function aplotStationGrid(dayRoots, stations, outDir)
% APLOTSTATIONGRID  Station columns x metric rows (PRDS/YEL3/SCOR style):
%   row 1: DTEC median per station + RMS annotation
%   row 2: cycle-slip counts + Total annotation
%   row 3: ROTI per-satellite colored dots
% Uses first available day per station (or pools all days if spanDays=true).
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(dayRoots{1}, '..', 'ArticleFigs');
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    nS = numel(stations);
    fig = figure('Name', 'Station grid', ...
        'Position', [100 50 380*nS + 120, 900], 'Visible', 'off');
    cmap = lines(32);
    for s = 1:nS
        % pool all days for this station
        Dall = []; Rall = []; Mall = [];
        for d = 1:numel(dayRoots)
            S = aloadTEC(dayRoots{d}, stations{s});
            if isempty(S) || ~isfield(S, 'TEC'), continue; end
            if isfield(S, 'obs')
                [~, MW, ~, ~] = agf_mw(S.obs);
                Mall = cat(1, Mall, MW); %#ok<AGROW>
            end
            Dall = cat(1, Dall, S.TEC.vertical); %#ok<AGROW>
            Rall = cat(1, Rall, S.ROTI); %#ok<AGROW>
        end
        if isempty(Dall)
            continue;
        end
        tH = (0:size(Dall, 1)-1) / 3600;
        % row 1: DTEC median + RMS
        ax = subplot(3, nS, s); hold on; grid on;
        dm = nanmedian(Dall, 2);
        dd = adtec(dm, 3600);
        plot(tH, dd, 'b-', 'LineWidth', 1);
        ok = isfinite(dd);
        rmsv = sqrt(nanmean(dd(ok).^2));
        text(0.05, 0.85, sprintf('RMS: %.3f m', rmsv), ...
            'Units', 'normalized', 'FontSize', 10, 'Color', 'b', 'FontWeight', 'bold');
        title(stations{s});
        ylabel('DTEC (TECU)');
        xlim([tH(1), tH(end)]);
        % row 2: cycle slips
        ax = subplot(3, nS, nS + s); hold on; grid on;
        nB = max(1, floor(size(Mall, 1) / 300));
        slips = zeros(nB, 1);
        for p = 1:size(Mall, 2)
            m = Mall(:, p);
            iv = find(isfinite(m));
            if numel(iv) < 2, continue; end
            dmw = abs(m(iv(2:end)) - m(iv(1:end-1)));
            idx = iv(2:end);
            idx = idx(dmw > 0.75);
            b = min(nB, floor((idx-1)/300)+1);
            for k = 1:numel(b), slips(b(k)) = slips(b(k)) + 1; end
        end
        bar((0:nB-1)*5/60, slips, 1, 'FaceColor', [0 0 0.8], 'EdgeColor', 'none');
        text(0.05, 0.85, sprintf('Total: %d', sum(slips)), ...
            'Units', 'normalized', 'FontSize', 10, 'Color', 'b', 'FontWeight', 'bold');
        ylabel('Num. of CS');
        xlim([tH(1), tH(end)]);
        % row 3: ROTI per sat
        ax = subplot(3, nS, 2*nS + s); hold on; grid on; %#ok<NASGU>
        keep = 1:20:size(Rall, 1);
        for p = 1:size(Rall, 2)
            v = Rall(keep, p);
            ok = ~isnan(v);
            if any(ok)
                plot((keep(ok)-1)/3600, v(ok), '.', 'Color', cmap(p, :), 'MarkerSize', 3);
            end
        end
        xlabel('UT (hour)');
        ylabel('ROTI (TECU/min)');
        xlim([tH(1), tH(end)]);
    end
    out = fullfile(outDir, 'aplotStationGrid.png');
    saveas(gcf, out);
    close;
    fprintf('aplotStationGrid saved: %s\n', out);
end
