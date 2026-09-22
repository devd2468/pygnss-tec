function [obs,nav,doy,Year] = readrinex304(r_o_name,r_n_name,rinex_path)
% Read observation and navigation RINEX file V. 3.02 - 3.04
% Inputs  : 
%          r_o_name - Observation RINEX's name
%          r_n_name - Navigation RINEX's name
% Outputs :
%          ----------- obs - observation ---------------
%          obs.date    - date, month, year
%          obs.epoch   - epoch of day (second of day)
%          obs.type    - type of pseudorange
%          obs.station - station name
%          obs.data    - data (pseudorange, SNR, etc.)
%          obs.index   - number of satellites
%          obs.rcvpos  - receiver position
%          ----------- nav - navigation ----------------
%          nav.eph     - ephemaris data of satellites
%          nav.index   - number of satellites
%          nav.ionprm  - Ionospheric coefficients
% You can change IGS station to download nav file
% Look at : ftp://anonymous:anonymous@cddis.gsfc.nasa.gov/pub/gps/data/daily/'year'/'doy'/'*n'/
%
% ****** If there are many lost data, the RINEX files cannot be read. ******

n_name = 'cnmr'; %amu2, ankr, cnmr, chpi, cusv
USE_MEX = false; % false = pure-MATLAB parsers (MEX binaries crash on R2026a)
current_path = [pwd '\'];
cd(rinex_path)
% Read RINEX Observation Files
if USE_MEX
    [Obs.Com,Obs.XYZr,Obs.Type,Obs.Date,Obs.Ep,Obs.Data,Obs.PRN] = ReadOBSrinex304(r_o_name);
else
    Obs = rinex3obs_gps(fullfile(rinex_path, r_o_name));
end
obs.date = round(Obs.Date.St(:)');
obs.epoch = Obs.Ep.G';
obs.type = Obs.Type.G;
obs.data = Obs.Data.G';
obs.index = Obs.PRN.G';
obs.rcvpos = Obs.XYZr;
obs.station = r_o_name(1:4);
% Multi-GNSS: keep GPS (.G) as primary (validated TEC pipeline) and, when the
% MEX reader provides other constellations, store them guarded (no crash if
% absent). Supported: G GPS, R GLONASS, C BeiDou, E Galileo, J QZSS.
obs.sys = struct();
try
    flds = {'R','E','C','J'};
    for k = 1:numel(flds)
        s = flds{k};
        if isfield(Obs.Type,s) && isfield(Obs.Ep,s) && isfield(Obs.Data,s) && isfield(Obs.PRN,s)
            obs.sys.(s).type  = Obs.Type.(s);
            obs.sys.(s).epoch = Obs.Ep.(s)';
            obs.sys.(s).data  = Obs.Data.(s)';
            obs.sys.(s).index = Obs.PRN.(s)';
        end
    end
catch
end
% Log constellations declared in the RINEX header (safe text parse, no MEX)
try
    obs.headerSys = detectRinexConstellations(fullfile(rinex_path, r_o_name));
    if ~isempty(obs.headerSys)
        fprintf('  Constellations in header: %s\n', strjoin(obs.headerSys, ','));
    end
catch
    obs.headerSys = {};
end

%     [obs.date, obs.epoch, obs.type, ~, ~, obs.station, obs.data, obs.index, obs.rcvpos,...
%         ~, ~, ~] = readrinexobs(r_o_name);

% Calculate Date Of Year (DOY)
disp('The data date:');
disp(datetime(obs.date(1:3)));
D1  = obs.date(1:3); % obs.date
D2  = D1; D2(:,2:3) = 0;
ydoy = cat(2, D1(:,1), datenum(D1) - datenum(D2));
Year = num2str(ydoy(1,1));
doy  = num2str(ydoy(1,2),'%.3d');
    
% Check NAV file \ Get online ephemeris
% NOTE (multi-GNSS/streamlined repo): if caller passes an explicit nav file
% that exists, use it directly and skip the legacy filename-pattern search
% (the old '*XXXX.XXn' pattern does not match long IGS names like
%  BHPL00IND_R_20251210000_01D_30S_MO.rnx) and skip auto-download.
useExplicitNav = ~isempty(r_n_name) && (exist(r_n_name,'file') || exist(fullfile(rinex_path,r_n_name),'file'));
if useExplicitNav
    if ~exist(r_n_name,'file')
        r_n_name = fullfile(rinex_path,r_n_name);
    end
    n = struct('name', r_n_name);
else
n = dir([rinex_path '*' r_o_name(end-7:end-1) 'n']);
end

if isempty(n) % Download Navigation file
    try
        nav_filename = [n_name doy '0.' Year(3:4) 'n.Z'];
        disp(['Download NAV RINEX file ' nav_filename])
        % download via built-in websave + gunzip (portable: no curl/gzip exes)
        websave(nav_filename, ['ftp://anonymous:anonymous@cddis.gsfc.nasa.gov/pub/gps/data/daily/' Year '/' doy '/' Year(3:4) 'n/' nav_filename]);
        gunzip(nav_filename);
        n = dir([rinex_path nav_filename(1:end-2)]);
        r_n_name = n.name;
        disp(['Nav file: ' r_n_name ' is downloaded'])
    catch
        cd(current_path)
        error(['error to download Nav file: ' n_name '. Please edit new Nav name in **readrinex304**'])
    end
end

% Read RINEX Navigation File
if ~USE_MEX
    navFull = r_n_name;
    if ~exist(navFull, 'file')
        navFull = fullfile(rinex_path, r_n_name);
    end
    Nav = rinex3nav_gps(navFull);
else
% --- Check D-, D+ --> E-, E+
Temp ='Temp.n';
Fnav = fileread(r_n_name); % read file (string)
if contains(Fnav,'D-')||contains(Fnav,'D+')
    nFnav = strrep(Fnav,'\','/');
    nFnav = strrep(nFnav,'D-','E-');
    nFnav = strrep(nFnav,'D+','E+');
    nFnav = strsplit(nFnav,'\n');
    [~,Col] = size(nFnav);
    nFid = fopen(Temp,'wt');
    [~,n] = size(nFnav{1});
    if n>80
        for i = 1:Col
            fprintf(nFid,[nFnav{i}(1:end-1) '\n']);
        end
    else
        for i = 1:Col
            fprintf(nFid,[nFnav{i} '\n']);
        end
    end
    fclose(nFid);
    clear nFnav
end
clear Fnav

[Nav.Com,Nav.Ion,Nav.dTime,Nav.Leap,Nav.Eph,Nav.PRN] = ReadNAVrinex304(Temp);
end
nav.eph = Nav.Eph.G';
nav.index = Nav.PRN.G';
nav.ionprm = Nav.Ion;
% Guarded multi-GNSS nav store (orbit models per system are future work)
try
    nav.sys = struct();
    for k = {'R','E','C','J'}
        s = k{1};
        if isfield(Nav.Eph,s) && isfield(Nav.PRN,s)
            nav.sys.(s).eph = Nav.Eph.(s)';
            nav.sys.(s).index = Nav.PRN.(s)';
        end
    end
catch
end
%     [~, ~, nav.eph, nav.index, nav.ionprm,...
%         ~, ~] = readrinexnav(r_n_name);

cd(current_path)
end

function sysList = detectRinexConstellations(obsFullPath)
% Safe header-only parse: return constellation letters from SYS / # / OBS TYPES
sysList = {};
try
    fid = fopen(obsFullPath, 'r');
    if fid < 0, return; end
    c = 0;
    while ~feof(fid) && c < 60
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        c = c + 1;
        if ~isempty(strfind(ln, 'SYS / # / OBS TYPES'))
            sysList{end+1} = strtrim(ln(1)); %#ok<AGROW>
        end
        if ~isempty(strfind(ln, 'END OF HEADER'))
            break;
        end
    end
    fclose(fid);
    % Keep only G,R,C,E,J (drop S/IRNSS per project scope), stable order
    keep = {'G','R','C','E','J'};
    sysList = intersect(keep, sysList, 'stable');
catch
    sysList = {};
end
end
