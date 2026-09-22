function cfg = config(varargin)
% CONFIG  Single configuration point for the multi-GNSS TEC pipeline.
% Works on any machine: all paths are absolute or anchored here; the only
% per-system edit a new user makes is the LOCAL PATHS block below (or pass
% name/value overrides, e.g. config('dataRoot', '/data/igs')).
%
% Fields:
%   dataRoot       folder with RINEX obs/nav (*_MO.rnx, *_MN.rnx, *.??o)
%   mgexRoot       folder with CODE MGEX SP3/CLK/ERP/OSB ('' = skip)
%   outRoot        results folder ('' = <dataRoot>/Results)
%   dateRange      {'YYYY-MM-DD', ...} day folders under dataRoot's parent
%                  ('' = process dataRoot itself as a single day)
%   stations       {} = all found, or {'BHPL','LCK4',...}
%   constellations {'G','R','C','E','J'} subset to attempt
%   elevMask       degrees (default 30)
%   ionoHeightKm   thin-shell height (default 350)
%   rotiMinor/rotiModerate/rotiSevere  TECU/min thresholds
%   minEventMin / maxGapMin            event detection (minutes)
%   allowDownload  false = never touch the network (default true-off safe)
%   usePrecise     false for now (broadcast ephemeris; SP3/CLK staged later)

%% ============ LOCAL PATHS (edit for your machine) ============
cfg.dataRoot = 'D:\IGS DATA\May 2025\2025-05-01';
cfg.mgexRoot = 'D:\IGS DATA\May 2025\MGEX';
cfg.outRoot  = '';   % '' = <dataRoot>/Results

%% ============ SELECTIONS ============
cfg.dateRange      = {};                 % e.g. {'2025-05-01','2025-05-05'}
cfg.stations       = {};                 % e.g. {'BHPL','LCK4'}
cfg.constellations = {'G', 'R', 'C', 'E', 'J'};
cfg.elevMask       = 30;
cfg.ionoHeightKm   = 350;
cfg.rotiMinor      = 0.2;
cfg.rotiModerate   = 0.5;
cfg.rotiSevere     = 1.0;
% S4 amplitude scintillation thresholds (standard literature values)
% Conker et al. (2003); Van Dierendonck et al. (1993)
cfg.s4Weak         = 0.15;   % background / low activity
cfg.s4Moderate     = 0.30;   % moderate amplitude scintillation
cfg.s4Strong       = 0.50;   % strong scintillation (may cause loss-of-lock)
cfg.minEventMin    = 15;
cfg.maxGapMin      = 5;
cfg.allowDownload  = false;
cfg.usePrecise     = false;

%% ============ Novel Research Analysis (Q1 GPS Solutions) ============
cfg.runNovelAnalysis = false;  % set true to execute all 5 novel modules
%  Solar-weather for novelAnalysis_SolarMax (1-5 May 2025 defaults):
cfg.F107 = [155, 158, 162, 170, 165];  % daily F10.7 index (sfu)
cfg.Kp   = [1.5,  1.8,  1.2,  2.3, 1.0];  % daily Kp index

%% ============ overrides: config('dataRoot', X, 'stations', {...}) ============
for k = 1:2:numel(varargin)
    if k+1 <= numel(varargin) && isfield(cfg, varargin{k})
        cfg.(varargin{k}) = varargin{k+1};
    end
end

%% ============ derived + validated ============
if isempty(cfg.outRoot)
    cfg.outRoot = fullfile(cfg.dataRoot, 'Results');
end
if ~exist(cfg.dataRoot, 'dir')
    error('config:dataRoot', 'dataRoot not found: %s', cfg.dataRoot);
end
if ~isempty(cfg.mgexRoot) && ~exist(cfg.mgexRoot, 'dir')
    warning('config:mgex', 'mgexRoot not found, continuing without MGEX: %s', cfg.mgexRoot);
    cfg.mgexRoot = '';
end
validSys = {'G', 'R', 'C', 'E', 'J'};
cfg.constellations = intersect(cfg.constellations, validSys, 'stable');
if isempty(cfg.constellations)
    error('config:constellations', 'No valid constellations selected.');
end
end
