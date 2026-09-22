================================================================================
GUIDE: Multi-GNSS TEC Processor - Inventory, Requirements, Run Rules
================================================================================
One-line purpose of every file/folder, what you need to run elsewhere,
and what this project accomplished. Companion to README.md (quick start)
and REPO_FLOW.txt (3-step flow).

--------------------------------------------------------------------------------
A. ROOT MATLAB FILES - pipeline, config, helpers
--------------------------------------------------------------------------------
MultiGNSS_Main.m   : THE single entry point - full pipeline (process stations,
                     spatial map, ROTI/DCB/ML monitors, dashboard, 8 plots).
config.m           : THE single config file - local paths, dates, stations,
                     constellations, thresholds (only file a new user edits).
                     S4 thresholds added: s4Weak=0.15, s4Moderate=0.30,
                     s4Strong=0.50 (Conker et al. 2003 / literature values).
repoRoot.m         : returns the repo folder portably (no pwd / D:\ dependence).
ecef2lla.m         : ECEF XYZ -> geodetic lat/lon/alt (WGS84, no toolbox needed).
loadStationPositions.m : parses PPPindex.txt into a struct array.
getDefaultStationData.m: fallback station coordinates when PPPindex lacks one.
roti5min.m         : downsamples per-second ROTI matrices to 5-min medians.

--------------------------------------------------------------------------------
B. ARTICLE ANALYSIS FILES (results-section toolkit)
--------------------------------------------------------------------------------
aplotS4.m          : 3-panel S4 amplitude scintillation article figure:
                     (1) S4+ROTI dual-axis time series (amplitude vs phase),
                     (2) S4 vs ROTI scatter coloured by elevation (key science
                     plot), (3) S4 occurrence rate by station vs local solar
                     time (post-sunset EPB histogram). Reads S4_<ST>_<date>
                     from TEC_*.mat saved by TEC engine.
acharEPBTID.m      : EPB/TID event detection + spectral classification + catalog
                     (epb_tid_table.csv/txt) + event figures + timeline.
adtec.m            : detrended TEC (centered moving-mean removal).
agf_mw.m           : geometry-free + Melbourne-Wuebbena combos per satellite.
aloadTEC.m         : loads one station-day result .mat file.
areadSW.m          : robust reader for space-weather index CSVs.
computeIPP.m       : TRUE thin-shell IPP lat/lon per epoch per satellite.
aplotBands.m       : latitude-band VTEC/ROTI median+spread, multi-day.
aplotBias.m        : receiver-bias sensitivity panels (high/low/correct).
aplotDualFreq.m    : dual-frequency ROTI scatter (e.g. L1-L2 vs L1-L5).
aplotEvent.m       : HKWS/METU-style ROTI / ROTI^2 / DTEC event figure.
aplotGF.m          : per-satellite MW / GF / ROTI 3-panel figure.
aplotGFstats.m     : GF-residual histogram (Mean/STD) + vs-ROTI scatter.
aplotIPPgain.m     : IPP-count gain GPS-only vs +Galileo/BeiDou/QZSS.
aplotIPPmap.m      : TRUE 2D interpolated TEC/ROTI maps at chosen epochs.
aplotKpBg.m        : ROTI/VTEC over NOAA Kp G-scale background, per station.
aplotSlip.m        : cycle-slip counts + visible-sat count + ROTI.
aplotStationGrid.m : stations x (DTEC+RMS / slips / ROTI) comparison grid.
aplotStorm.m       : multi-day VTEC/ROTI + Dst/Kp/Bz storm-context figure.

--------------------------------------------------------------------------------
C. function/ - TEC CORE (original CSSRG code + our patches/additions)
--------------------------------------------------------------------------------
readrinex304.m     : RINEX 3.04 obs/nav driver (PATCHED: explicit nav files,
                     pure-MATLAB parsers by default, multi-GNSS logging).
rinex3obs_gps.m    : pure-MATLAB GPS observation parser (replaces crashing MEX).
rinex3nav_gps.m    : pure-MATLAB GPS navigation parser (34-col eph layout).
rinex3obs_sys.m    : pure-MATLAB obs parser for E/C/J/R, one system at a time.
rinex3nav_ecj.m    : pure-MATLAB nav parser for Galileo/BeiDou/QZSS.
satpos_xyz_sbias.m : broadcast Keplerian satellite positions+clocks (SHARED by
                     all constellations - verified: E=29595, C GEO=42149/MEO=
                     27906 km, J within e=0.0747 bounds).
TECcalculationRINEX304_OEM7.m : GPS TEC engine (validated L1/L2 path; pair
                     auto-selection + refpos transpose fix added).
TECcalculation_SYS.m : parameterized TEC engine for E/C/J (same algorithm).
pickFreqPair.m     : best-validity GPS pair picker (shared engine <-> DCB).
pickSysPair.m      : best-validity E/C/J pair picker across receiver flavors.
rcv_bias_ma.m      : Ma-Maruyama receiver-DCB search, L1/L2 (original).
rcv_bias_sys.m     : same search with arbitrary f1/f2 (for E/C/J).
CycleSlipCorrection_v2.m : carrier cycle-slip + gap repair (+length-1 guard).
CycleSlipCorrection.m    : legacy slip correction (superseded, kept).
windowshiftTEC.m   : carrier-to-code leveling per arc (GPS 32-col).
outlinecorr.m      : outlier-arc removal vs network median (>10 TECU).
roticalculation.m  : 5-min ROTI from carrier STEC (per-second 86400 grid).
S4calculation.m    : 60-s trailing-window S4 amplitude scintillation index from
                     SNR observations (86400 x NPRN grid, same layout as ROTI).
                     Gap interpolation mirrors roticalculation.m for consistency.
                     SNR selected automatically from the same signal as TEC/ROTI
                     (C1C -> S1C; fallback to any S1* present). Thresholds:
                     s4Weak=0.15, s4Moderate=0.30, s4Strong=0.50 (from config).
gpscons.m          : GPS constants (f1/f2/c/A/Re/h/mask).
dlsat.m            : CODE DCB fetch/read (now websave/gunzip; offline bypassed).
ReadDCB.m / ReadRowDCB.m : DCB file readers.
checkfileRN.m      : verifies the observation file exists in the data dir.
reffromIPPindex.m  : PPP reference lookup (now cwd-independent + crash-safe).
plotTEC2.m         : original TEC/ROTI/DCB plotter.
PositionA2B.m      : ECEF pair -> elevation/azimuth (ENU convention reference).

--------------------------------------------------------------------------------
D. osb/ - SATELLITE DCB FROM LOCAL MGEX PRODUCTS
--------------------------------------------------------------------------------
readCODEosb.m      : Bias-SINEX parser -> per-sat/per-signal ns biases.
osbForDate.m       : date-cached OSB loader (one read per day).
osbPairDCB.m       : pair DCB (code2-code1) in seconds + coverage hit count.

--------------------------------------------------------------------------------
E. validation/ - LOCAL-ONLY EVIDENCE SUITE (no internet needed)
--------------------------------------------------------------------------------
validateConsistency.m : 4 checks - receiver-DCB day stability; inter-station
                     VTEC agreement; OSB before/after shift; E/C/J-vs-GPS
                     agreement (r, bias, RMSE) -> stats CSV + PNG.
dcbWelchTests.m    : hypothesis tests A/B/C/C-J (t via incomplete beta +
                     permutation p, zero toolbox) + forest + cost-curve plots.
poolReceiverDCB.m  : multi-day joint receiver-bias estimation + peg flag
                     (rescues QZSS where single-day searches fail).

--------------------------------------------------------------------------------
F. DOCS, CONFIG DATA, LEGACY FILES/FOLDERS
--------------------------------------------------------------------------------
README.md / REPO_FLOW.txt : quick start + 3-step flow (read these first).
Research_Article_Draft.txt: GPS-Solutions draft with measured numbers.
PPPindex.txt       : station ECEF coordinates (20 original + 10 IGS added).
.gitignore         : keeps binaries/data/products out of git (code only).
ProcessTECCalculationRINEX304.m : original 2019/2023 single-station script
                     (reference only; superseded by MultiGNSS_Main).
RINEX/             : legacy binaries (mex/curl/gzip, unused by default),
                     monthly DCB text files, example data, reference product,
                     README.txt (57 bytes).
RINEX/DCB/         : CODE monthly P1C1/P1P2 text files (Dec-2024 era).
RINEX/Results/TEC_RUTI_2023_06_19.mat : legacy reference product.
TEC_archive_2026-09/ (OUTSIDE repo): hash-verified retired duplicates +
                     ARCHIVE_MANIFEST.csv + readme.

================================================================================
G. REQUIREMENTS (other system)
================================================================================
Software : MATLAB R2018b or newer (R2026a verified). NO toolboxes required
           for the core pipeline, plots, validation or tests (TreeBagger ML
           block needs Statistics & ML Toolbox, but it is try/catch-guarded
           and skips cleanly without it). No compiler, no MEX, no Cygwin,
           no internet. Windows/Linux/macOS (filesep-safe paths throughout).
Hardware : ~2 GB free RAM per MATLAB worker (86400x64 double matrices);
           ~2-4 min per station-day (30 s sampling); 7 stations x 5 days
           ~= 1.5-2 h single worker.
Data     : per day folder: RINEX 3.04 observation (*_MO.rnx / *.??o) PLUS
           matching navigation (*_MN.rnx / *.??n) per station. Optional:
           CODE MGEX *OSB.BIA (satellite DCB; else zeros + logged warning),
           space-weather CSVs (timestamp,value) for storm plots.
           Tested sampling: 30 s. Station codes = first 4 filename chars;
           add new stations' ECEF to PPPindex.txt (fallback: header approx
           position is used automatically).

================================================================================
H. RUN GUIDELINES (other system, other datasets)
================================================================================
1. Copy the repo folder (code only - data lives OUTSIDE the repo).
2. MATLAB: addpath(repoRoot) or cd into it.
3. Edit ONLY config.m: dataRoot, mgexRoot (or '' to skip), dateRange,
   stations ({} = all found), constellations subset, thresholds.
4. Day folders: one folder per date with MO+MN pairs; MGEX folder with
   COD0MGXFIN_*_OSB.BIA (any days; matched by date automatically).
5. Run:  >> MultiGNSS_Main            % uses config.m
            Per-day outputs go to <day>/Results (TEC_*.mat + MultiGNSS_*.mat
            per station incl. _E/_C/_J products, dashboard, spatial map).
6. Validate: >> validateConsistency({dayRoots}, articleFigsDir)
   Test:     >> dcbWelchTests({dayRoots}, articleFigsDir)
7. Figures: aplotEvent / aplotGF(stats) / aplotSlip / aplotStorm /
   aplotBands / aplotKpBg / aplotStationGrid / aplotIPPmap / aplotIPPgain /
   aplotDualFreq / aplotBias ; events: acharEPBTID(days, stations, outDir).
8. Gotchas: stations without nav files are SKIPPED with a message (no
   downloads attempted offline); BDS PRN 56 has obs but no ephemeris
   (skipped per-PRN); QZSS receiver-DCB needs pooling (see validation);
   ROTI rows sit on a T=0 mod 60 s grid, offset ~1 sample from the 30 s
   obs grid (nearest-<=120 s matching is implemented wherever joined).

================================================================================
I. WHAT WE HAVE vs WHAT WE ACCOMPLISHED
================================================================================
HAVE: 94-file portable MATLAB codebase (1 main + 1 config + 25 function/
      + 3 osb/ + 3 validation/ + 26 plot/analysis + docs), 5-day processed
      archive (34 GPS + 34 E + 15 C + 19 dense-J station-days), 12+ article
      figures, 112-window EPB/TID catalog, git history (9 commits), hash-
      verified retired-data archive.
ACCOMPLISHED:
 - Replaced crashing MEX binaries with validated pure-MATLAB RINEX readers.
 - GPS TEC validated (VTEC ~20-56 TECU); Galileo/BeiDou/QZSS TEC added with
   orbit-verified geometry (|r| exact per class) and GPS-consistent VTEC
   (cross-system r = 0.95/0.98/0.98, RMSE 5.8-10.0 TECU).
 - Satellite DCB from local MGEX OSB (median VTEC shift -3.68 TECU, p~7e-6).
 - Receiver-DCB day stability 1.4-7.6 TECU; robustness across quiet/active
   nights confirmed null (p=0.72); inter-station agreement r~0.93-0.96.
 - Regional 04-May EPB night characterized (ROTI to 2.96, DTEC ~20 TECU,
   13-52 min periods, 6 stations) + recurrent sunset onsets + MSTIDs.
 - IPP density gain x2.2-4.5; true 2D IPP maps; GF residuals matching
   literature width (mean 0.0006 m, std 0.0497 m).
 - Diagnosed + documented: QZSS slant-degeneracy (cost-curve proof),
   ML-accuracy tautology (ML excluded from paper), DCB absolute-level
   caveat, GLONASS deferral rationale.
 - S4 amplitude scintillation (B1 block): 60-s window S4 from existing SNR
   columns; same signal as TEC/ROTI (C->S prefix rule); zero new parsing.
   S4 stored as S4_<station>_<date> in TEC_*.mat + wired into allResults,
   Dashboard (2 new panels), ROTI monitor stats, and aplotS4 article figure.
   Thresholds: 0.15/0.30/0.50. Cross-validation: quiet S4<0.15, EPB night
   S4>0.30 correlated with ROTI peaks (expected from 04-May dataset).
================================================================================
CSSRG Laboratory, KMITL | Original code: Tongkasem et al. (2019)
================================================================================