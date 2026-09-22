function [time_h, med5] = roti5min(ROTI)
% ROTI5MIN  Downsample ROTI matrix to 5-minute medians for display.
% ROTI from roticalculation.m is 86400x32 (per-second rows); legacy files
% may be 288x32 (5-min rows). Returns time in hours + median series.
    nR = size(ROTI, 1);
    if nR > 1000
        nb = floor(nR / 300);
        med5 = nan(nb, 1);
        for b = 1:nb
            blk = ROTI((b-1)*300+1:min(b*300, nR), :);
            med5(b) = nanmedian(blk(:));
        end
        time_h = (0:nb-1) * 5 / 60;
    else
        med5 = nanmedian(ROTI, 2);
        time_h = (0:nR-1) * 5 / 60;
    end
end
