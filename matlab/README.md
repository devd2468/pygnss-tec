================================================================================
                    CSSRG TEC Research Project
                    Multi-GNSS TEC Processor - Streamlined Repo
================================================================================

Single main file: MultiGNSS_Main.m
  GPS (G) + Galileo (E) + BeiDou (C) + QZSS (J) TEC  [R/GLONASS: detected only]
  + broadcast nav + MGEX OSB satellite DCB + spatial maps + ROTI monitor
  + DCB stability + article plots (event/GF/storm/bands/Kp/slips/grid)

================================================================================
REPO STRUCTURE
================================================================================

RUN THIS:
  MultiGNSS_Main.m          <- single entry point (all functions inside)

HELPERS (called by main / plots):
  ecef2lla.m, loadStationPositions.m, getDefaultStationData.m
  plotMagneticEquator.m, plotMappingRange.m, plotMapWithIPP.m,
  plotNortheastReceivers.m, plotROTIdtrend.m, plotROTIinDifferentsSites.m,
  plotROTImaps.m, plotScatterAndInterpolation.m

ORIGINAL CODEBASE (untouched, in function/):
  readrinex304.m  (patched: explicit-nav + guarded multi-GNSS fields)
  TECcalculationRINEX304_OEM7.m, roticalculation.m, gpscons.m,
  satpos_xyz_sbias.m, CycleSlipCorrection_v2.m, windowshiftTEC.m,
  rcv_bias_ma.m, outlinecorr.m, plotTEC2.m, dlsat.m, ReadDCB.m,
  ReadRowDCB.m, checkfileRN.m, reffromIPPindex.m, PositionA2B.m, ...

DATA / DOCS:
  PPPindex.txt, README.md, REPO_FLOW.txt, Research_Article_Draft.txt
  RINEX/  (mex + curl/gzip + example data + Results)

================================================================================
QUICK START
================================================================================

1. EDIT config.m (one file: local paths, date range, stations,
   constellations, thresholds), or pass overrides:
     >> MultiGNSS_Main
     >> MultiGNSS_Main('/data/igs/2025-05-01')
     >> MultiGNSS_Main('dataRoot', X, 'stations', {'BHPL'})
   Works on any OS / any folder - no hardcoded paths in code.

2. VIEW results in <DATA>\Results\ (dashboard PNG, spatial map, 8 plots,
   ROTI_events.mat, DCB analysis, Phase3 model).

3. VALIDATE (local-only, no internet):
     >> validateConsistency({dayRoot1, dayRoot2, ...}, articleFigsDir)
   Checks: receiver-DCB day stability; inter-station VTEC agreement;
   OSB before/after shift; Galileo/BeiDou/QZSS-vs-GPS agreement
   (r, bias, RMSE). Results: validation_stats.csv + PNG.

================================================================================
NOTES
================================================================================
- Observation filter: only *_MO.rnx / *.??o are processed; *MN/*EN/*GN/*RN
  navigation files are NEVER fed to the observation MEX reader (this was the
  crash cause).
- TEC engines: validated GPS L1/L2 path + Galileo E1/E5a(E5b), BeiDou
  B1I/B2I (B1C/B2a fallback), QZSS L1/L2(L5) via parameterized
  TECcalculation_SYS (same algorithm family). GLONASS deferred (FDMA +
  orbit integrator out of scope).
- Verified orbits: |r| E=29595 km, C GEO=42149/MEO=27906 km, J matches
  e=0.0747 Keplerian bounds. IPP gain x2.2-4.5 vs GPS-only.
- Satellite DCB from local MGEX OSB (P1-style: OSB_C1W-OSB_Cx); E5b/B2I
  pairs lack OSB coverage -> zeros (logged per station).
- Pure-MATLAB RINEX readers (no MEX binaries); downloads use built-in
  websave/gunzip (no curl.exe, gzip.exe, or Cygwin needed) - portable.
================================================================================
CSSRG Laboratory, KMITL | Original: Napat Tongkasem et al. (2019)
================================================================================