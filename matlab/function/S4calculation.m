function S4_sec = S4calculation(SNR_mat, PRNall)
%{
 ========================================
           S4 Amplitude Scintillation Index (60-s trailing window)
 ========================================

 Description
 -----------
 Computes the S4 amplitude scintillation index from GNSS Signal-to-Noise
 Ratio (SNR) observations. S4 is defined as:

     S4 = std(SNR) / mean(SNR)  over a 60-second trailing window

 This mirrors the structure and gap-interpolation logic of roticalculation.m
 so that S4 and ROTI are computed on the same per-second 86400-row grid
 and can be directly compared (amplitude vs. phase scintillation).

 SNR is taken from the obs.data columns labelled S1C / S2W / S5X (or
 whichever S-type corresponds to the selected code pair). The readers leave
 SNR values untouched (zero-mapping applies only to C*/L* columns), so
 raw dB-Hz values are used here.

 Inputs
 ------
   SNR_mat  - (86400 x NPRN) SNR matrix, one column per PRN; NaN where
              no observation. Values in dB-Hz (receiver-reported S* obs).
   PRNall   - row vector of PRN indices (1-based) to process.

 Output
 ------
   S4_sec   - (86400 x NPRN) S4 index on the same per-second grid.
              NaN where fewer than 3 valid SNR samples exist in the window.

 Thresholds (from config.m)
 --------------------------
   s4Weak     = 0.15  (background / low activity)
   s4Moderate = 0.30  (moderate amplitude scintillation)
   s4Strong   = 0.50  (strong; comparable to S4 > 0.5 in literature)

 References
 ----------
   Conker, R.S., et al. (2003). Modeling the effects of ionospheric
     scintillation on GPS/SBAS availability. Radio Sci., 38(1), 1001.
   Fremouw, E.J., et al. (1978). Early results from the DNA Wideband
     satellite experiment. Radio Sci., 13, 167-187.
   Van Dierendonck, A.J., et al. (1993). Theory and performance of
     narrow correlator spacing in a GPS receiver. NAVIGATION, 39, 265-283.

 CSSRG Laboratory, KMITL | S4 implementation: B1 (S4) block, Sep 2026
%}

S4_sec = nan(size(SNR_mat));   % pre-fill with NaN (86400 x NPRN)

% Window: 60-second trailing window (60 per-second samples)
WIN = 60;

for PRN = PRNall
    SNR_col = SNR_mat(:, PRN);

    % ---- Gap interpolation (identical to roticalculation.m lines 17-28) ----
    % Short gaps (1 < gap <= 300 s) are bridged by linear interpolation
    % to avoid artificial S4 spikes at arc boundaries. This keeps S4 on the
    % same arc-continuity assumption as ROTI.
    ST = find(~isnan(SNR_col));
    if isempty(ST)
        continue;          % no data for this PRN
    end

    flag = find((diff(ST)) > 1 & (diff(ST)) <= 300);   % short-gap epochs
    if ~isempty(flag)
        for d = 1:length(flag)
            x  = [1, length(SNR_col(ST(flag(d)):ST(flag(d)+1)))];
            v  = [SNR_col(ST(flag(d))), SNR_col(ST(flag(d)+1))];
            xq = 1:length(SNR_col(ST(flag(d)):ST(flag(d)+1)));
            SNR_M = interp1(x, v, xq, 'linear', 'extrap');
            SNR_col(ST(flag(d)):ST(flag(d)+1)) = SNR_M;
        end
    end

    % ---- S4 over 60-s trailing window ----
    for T = WIN:length(SNR_col)
        if isnan(SNR_col(T))
            S4_sec(T, PRN) = NaN;
            continue;
        end

        % Extract window; keep only finite values
        window = SNR_col(T-WIN+1:T);
        valid  = window(isfinite(window));

        if length(valid) < 3
            % Fewer than 3 valid samples -> unreliable estimate, leave NaN
            continue;
        end

        mu_snr = nanmean(valid);
        if mu_snr < 1e-6
            % Near-zero mean SNR: guard against divide-by-zero
            continue;
        end

        % Standard S4 definition: sigma(I)/mean(I)
        % Using SNR directly as a proxy for signal intensity I.
        % For dB-Hz data this is an approximation; for linear-power SNR it is
        % exact. The ratio is dimensionless in either case, and relative
        % comparisons between epochs/satellites remain valid.
        S4_sec(T, PRN) = nanstd(valid) / mu_snr;
    end
end
end
