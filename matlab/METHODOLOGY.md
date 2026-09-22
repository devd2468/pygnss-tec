================================================================================
METHODOLOGY & WORKFLOW: Multi-GNSS Ionospheric TEC/ROTI Analysis
================================================================================
How the analysis is done, step by step, and everything that is assumed,
thresholded, or decided along the way. Companion to GUIDE.md (inventory)
and REPO_FLOW.txt (3-step run guide).

--------------------------------------------------------------------------------
1. PHYSICAL MODEL
--------------------------------------------------------------------------------
Dual-frequency geometry-free combination per satellite, per epoch:
    STECp = k * (P2 - P1)          % from code pseudoranges (m), noisy
    STECl = k * (L1 - L2)          % from carrier phase x wavelength (m), precise
                                   %   but ambiguous (unknown integer offset)
with k = f1^2*f2^2 / (A*(f1^2-f2^2)*1e16)  [TECU per meter], A = 40.3,
c = 299792458 m/s. k is recomputed per frequency pair:
    GPS L1/L2 (1575.42/1227.60 MHz): k = 9.5196
    Galileo E1/E5a, BeiDou B1C/B2a, QZSS L1/L5: k ~= 7.7637
    Galileo E1/E5b: k ~= 8.7575 | BeiDou B1I/B2I (1561.098/1207.14): k ~= 8.9932
Wavelengths lam = c/f convert cycles to meters. Only L1/L2-class pairs are
used (same k-form); triple-frequency combinations are not implemented.

--------------------------------------------------------------------------------
2. STEP-BY-STEP PIPELINE (per station-day)
--------------------------------------------------------------------------------
S0 INTAKE. Day folder holds RINEX 3.04 observation (*_MO.rnx / *.??o) plus
   matching navigation (*_MN.rnx / *.??n). Files whose names match
   navigation patterns are NEVER fed to the observation reader (this was the
   historic crash cause). Stations without a local nav file are SKIPPED with
   a message; no downloads are attempted offline (config.allowDownload=false).

S1 PARSE (pure MATLAB, no binaries). Observation: header SYS types per
   constellation (multi-line continuations handled), APPROX position,
   TIME OF FIRST OBS; body epochs (flag>1 skipped), 16-char fields,
   blank->NaN. Loss-of-lock written as 0.000 is mapped to NaN for C*/L*
   (0.0 code/phase is never physical; Doppler/SNR untouched). Navigation:
   8-line GPS-like Keplerian records (GPS/Galileo/BeiDou/QZSS), Fortran
   D-exponents, measured column slices (line1 Toc 4:23/data 24,43,62;
   lines 2-8 data 5,24,43,62). GLONASS 4-line records are parsed for
   detection only (orbits deferred).

S2 ORBITS. Unmodified broadcast-Kepler solver (satpos_xyz_sbias): iterative
   Kepler equation (1e-12 tolerance), J2-perturbation corrections, Earth-
   rotation (We = 7.2921151467e-5 rad/s, GM = 3.986004418e14),
   relativistic clock term. Galileo/BeiDou/QZSS ephemerides are decoded
   into the identical 34-column layout; their weeks/TOWs are already on the
   GPS SOW grid in these files (verified: BDS Toe=342000 at Wed 23:00).
   Verified evolving mean radius: Galileo ~29,595 km, BDS GEO ~42,149 km,
   BDS MEO ~27,906 km, QZSS within e=0.0747 bounds.

S3 SIGNAL SELECTION. Best-validity L1/L2-class pair per file (shared
   pickFreqPair/pickSysPair used identically by TEC engine and DCB builder):
   GPS C1C/L1C+C2W/L2W or C2X/L2X; Galileo E1/E5a (E1/E5b fallback);
   BeiDou B1I/B2I (B1C/B2a fallback); QZSS L1/L2 (L1/L5 fallback).

S4 GEOMETRY. Elevation from ECEF receiver-satellite vectors; mask < 30 deg
   (gpscons elev_mask). Reference position: PPPindex.txt hit, else the
   file's own APPROX position (row-vector guarded). Thin shell h = 350 km;
   slant factor sqrt(1-(Re*cos(el)/(Re+h))^2), Re = 6371.009 km.

S5 CARRIER PROCESSING. Cycle-slip/missing-data repair
   (CycleSlipCorrection_v2; slip threshold = std of STEC steps clamped to
   1-2 TECU; arcs < 2 samples pass through untouched), carrier-to-code
   leveling per continuous arc (windowshiftTEC; re-split after >=3 h gaps),
   zero-floor shift, outlier-arc rejection vs network median (>10 TECU,
   outlinecorr), tail pruning (final ~25 min, epochs 84900-86400, dropped
   per PRN if <1500 valid samples there).

S6 BIASES. Satellite DCB priority: local CODE MGEX OSB (per-signal ns
   biases, cached per date; P1Cx = OSB_C1W - OSB_Cx; missing sats = 0 and
   counted) -> legacy monthly DCB files if present -> zeros. Receiver DCB
   by Ma-Maruyama minimum-spread search (-30..+30 ns, 0.1 step, 6 refinement
   passes, 1:30 downsampling). Final: STEC -= satDCB; VTEC = (STEC-rcvDCB)*
   slant, zero-shifted to minimum.

S7 ROTI (Pi et al. 1997). Per-second grid: at each row T>=300, std of up to
   six 1-min carrier samples (T-299:60:T). NOTE (documented grid offset):
   values land on rows T = 0 mod 60 s, ~1 sample off the 30 s obs grid;
   all joins use nearest-valid <= 120 s. Thresholds 0.2/0.5/1.0 TECU/min.

S7b S4 AMPLITUDE SCINTILLATION (B1 block — companion to S7/ROTI).
   Input: SNR observations from obs.data columns (S* type). SNR is NOT
   zero-mapped by the reader (only C*/L* are zero->NaN; SNR kept raw).
   Signal selection: derive SNR type from the selected code pair by prefix
   substitution (C->S): tC1='C1C' -> snrType='S1C', etc. Falls back to any
   S1* column on the same frequency if the exact match is absent. This
   ensures S4 is on the same signal as TEC and ROTI (apples-to-apples
   amplitude-phase comparison essential for the science plot).
   Algorithm (S4calculation.m, mirrors roticalculation.m):
     1. Build 86400 x NPRN SNR matrix (same epoch grid as STECl), elevation
        mask applied (epochs outside 30 deg mask set to NaN).
     2. Gap interpolation: linear bridging of gaps 1 < Δt <= 300 s (same
        logic as roticalculation.m lines 17-28) to avoid arc-boundary spikes.
     3. At each second T >= 60: extract 60-sample trailing window; require
        >= 3 finite values; compute S4 = std(SNR_window)/mean(SNR_window).
   Output: 86400 x NPRN S4_sec matrix (NaN where insufficient data).
   Note on SNR units: receivers report S* in dB-Hz. Using dB-Hz directly is
   an approximation (true S4 uses linear power); the ratio std/mean is
   dimensionless and relative comparisons between epochs and satellites are
   internally consistent. For absolute S4 values, convert to linear first.
   Saved as S4_<station>_<date> (GPS) and S4_<tag>_<date> (E/C/J) inside
   TEC_*.mat alongside ROTI, and loaded into allResults.S4.
   Thresholds: s4Weak=0.15, s4Moderate=0.30, s4Strong=0.50 (config.m).
   Verification criteria:
     quiet periods  -> S4 < 0.15  (background noise floor)
     04-May EPB night -> S4 > 0.30 correlated with ROTI >= 0.5 peaks
     S4 finite only where SNR is valid (no synthetic values injected)

S8 PRODUCTS (per station-day, per system). 86400-row matrices (GPS 32 cols,
   E/C/J 64 cols): TEC.vertical/slant/withrcvbias/withbias/STECp/STECl,
   ROTI, S4 (amplitude scintillation, 60-s window, same grid as ROTI),
   DCB.sat/.rcv, prm.elevation, refpos. Saved twice: TEC_<ST>[_<SYS>]_
   <date>.mat (engine-native names, now includes S4_<ST>_<date>) +
   MultiGNSS_<ST>[_<SYS>]_<date>.mat (uniform TEC/ROTI/DCB/prm/S4 +
   obs/nav for GPS). allResults struct carries .S4 field for dashboard.

--------------------------------------------------------------------------------
3. POST-PROCESSING (multi-day, multi-station)
--------------------------------------------------------------------------------
P1 Spatial maps: station VTEC IDW/natural interpolation; true per-IPP maps
   (computeIPP recomputes satpos+azimuth from saved obs/nav, 350 km shell,
   elevation cross-check <0.1 deg) at chosen epochs (aplotIPPmap).
P2 MW/GF diagnostics: Melbourne-Wuebbena (cycles) + geometry-free (m) from
   saved obs; MW-jump>0.75 cyc slip proxy; GF residual stats + vs-ROTI
   scatter with envelope.
P3 EPB/TID catalog (acharEPBTID): contiguous per-PRN ROTI>=0.3 segments
   >=10 min (gaps <=5 min merged); per event: peak ROTI, DTEC (1-h detrend)
   peak-to-peak + depletion depth, FFT dominant period (10-120 min band),
   supporting-sat count. Rules: EPB = local 17.5-24 h IST AND (peak>=0.5
   or depth<=-5); TID = spectral 15-90 min peak with prominence>=2.
P4 Validation (local-only): receiver-DCB day stability; inter-station VTEC
   agreement (pairwise r of 5-min medians); OSB before/after shift;
   E/C/J-vs-GPS agreement (r/bias/RMSE, >=60 overlap bins); hypothesis
   tests A/B/C (t via incomplete beta + 10k permutation p, zero toolbox).

--------------------------------------------------------------------------------
4. WORKFLOW (execution order and data flow)
--------------------------------------------------------------------------------
config.m  -->  MultiGNSS_Main  -->  <day>/Results/*.mat
   (paths/dates/      (stations x systems:            (TEC/ROTI/DCB/prm,
    stations/systems/   parse -> orbit -> STEC ->       obs/nav for GPS)
    thresholds)         level -> DCB -> VTEC -> ROTI)
                                                        |
                          +-----------------------------+------------------+
                          |                             |                  |
                    validateConsistency          dcbWelchTests     aplot*/acharEPBTID
                    (4 checks -> CSV+PNG)        (A/B/C tests)     (article figures)
                                                        |
                                              ArticleFigs/*.png + tables
                                                        |
                                              Research_Article_Draft.txt

Single commands:
  >> MultiGNSS_Main                                     % config.m defaults
  >> MultiGNSS_Main('dataRoot', X, 'stations', {'BHPL'})
  >> validateConsistency({d1,d2,...}, figsDir)
  >> dcbWelchTests({d1,d2,...}, figsDir)
  >> evts = acharEPBTID({d1,...}, {'BHPL',...}, figsDir)

--------------------------------------------------------------------------------
5. KEY CONSIDERATIONS (what we account for, and known limits)
--------------------------------------------------------------------------------
- Absolute level: receiver DCB absorbs unmodeled common offsets; without
  external GIM the absolute VTEC carries several-TECU uncertainty (stated
  in paper; morphology/ROTI/DTEC unaffected).
- Ma-Maruyama nighttime-zero assumption is violated at solar-max EIA
  (receiver DCB biased high, precision 1.4-7.6 TECU/day retained).
- QZSS-class sparse/high-elevation geometries flatten the bias-search cost
  (boundary pegging) -> pooled multi-day estimation, descriptive-only stats.
- No per-epoch IPP azimuth in saved products -> recomputed on demand for maps.
- 30 s sampling assumed in lag/gradient code (30-row steps); ROTI grid note
  above; BDS PRN 56 has obs but no ephemeris (skipped per-PRN).
- Reproducibility: git history, hash-verified archive, per-run console log;
  randomness only in permutation tests (rng(42)) and plot decimation.
================================================================================
