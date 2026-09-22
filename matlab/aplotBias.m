function aplotBias(dayRoot, station4, outDir)
% APLOTBIAS  Receiver-bias sensitivity figure (Ma-Maruyama style).
%   (A) assumed receiver bias HIGH  (blue)
%   (B) assumed receiver bias LOW   (red)
%   (C) estimated receiver bias     (green)
% Each panel: per-satellite VTEC arcs + sigma_total annotation.
    if nargin < 3 || isempty(outDir), outDir = fullfile(dayRoot, 'Results'); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    S = aloadTEC(dayRoot, station4);
    if isempty(S) || ~isfield(S, 'TEC') || ~isfield(S, 'DCB')
        fprintf('aplotBias: no data for %s in %s\n', station4, dayRoot);
        return;
    end
    base = S.TEC.withrcvbias;   % STEC w/ receiver DCB, per PRN
    if ~isfield(S, 'prm') || ~isfield(S.prm, 'elevation')
        fprintf('aplotBias: no elevation for %s\n', station4);
        return;
    end
    Re = 6371.009e3; h = 350e3;
    sf = sqrt(1 - (Re * cosd(S.prm.elevation) ./ (Re + h)).^2);
    rcv = S.DCB.rcv;
    if ~isscalar(rcv) || isnan(rcv), rcv = 0; end

    delta = 15;  % TECU offset for high/low demonstration
    cases = struct('tag', {'(A) Assumed receiver bias is high', ...
                           '(B) Assumed receiver bias is low', ...
                           '(C) Estimated receiver bias is correct'}, ...
                   'off', {delta, -delta, 0}, ...
                   'col', {[0 0 0.7], [0.8 0 0], [0 0.6 0]});

    fig = figure('Name', sprintf('Bias %s', station4), ...
        'Position', [100 50 1000 1050], 'Visible', 'off');
    tH = (0:size(base, 1)-1) / 3600;
    for k = 1:3
        ax = subplot(3, 1, k); hold on; grid on;
        V = (base - (rcv + cases(k).off)) .* sf;
        for p = 1:size(V, 2)
            v = V(:, p);
            ok = ~isnan(v);
            if sum(ok) > 30
                plot(tH(ok), v(ok), '-', 'Color', cases(k).col, 'LineWidth', 0.8);
            end
        end
        sig = sigTotal(V);
        text(0.35, 0.75, sprintf('%s\n\\sigma_{Total} = %.2f, DCB_R = %.1f TECU', ...
            cases(k).tag, sig, rcv + cases(k).off), ...
            'Units', 'normalized', 'FontSize', 11, 'Color', cases(k).col);
        xlim([0 24]);
        ylabel('Vertical TEC (TECU)');
        if k == 3, xlabel('Time (UT)'); end
    end
    sgtitle(sprintf('%s  receiver-bias sensitivity', station4));
    out = fullfile(outDir, sprintf('aplotBias_%s.png', station4));
    saveas(gcf, out);
    close;
    fprintf('aplotBias saved: %s\n', out);
end

function s = sigTotal(V)
% total spread metric in the spirit of rcv_bias_ma (sum of per-sat std)
    s = 0;
    for p = 1:size(V, 2)
        v = V(:, p);
        v = v(~isnan(v));
        if numel(v) > 10
            s = s + std(v);
        end
    end
end
