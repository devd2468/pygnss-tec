function TECcalculation_SYS(obsS, navS, satBias_s, S_path, tag, F, PT)
% TECcalculation_SYS  TEC engine for one non-GPS constellation (E/C/J).
% Same validated algorithm as TECcalculationRINEX304_OEM7, parameterized:
% Inputs:
%        obsS  = struct(type 1xM cell, data NxM, index Nx1, epoch Nx1,
%                rcvpos 1x3, station string)
%        navS  = struct(eph Nrecx34 GPS-layout, index Nrecx1)
%        satBias_s = 64x1 satellite DCB (code2-code1) in seconds (0 unknown)
%        S_path    = results path (with trailing filesep)
%        tag       = save tag, e.g. 'BHPL_E_2025_05_01'
%        F     = struct(f1,f2,k,lam1,lam2) frequencies (Hz), k, wavelengths (m)
%        PT    = struct(c1,c2,l1,l2) observation code strings
% Saves: TEC_<tag>.mat with TEC_SYS_<date>, DCB_SYS_<date>, ROTI_SYS_<date>,
%        prm_SYS_<date>, refpos. Matrices are 86400x64 indexed by PRN.

% Ref position (same rule as GPS engine)
refpos = reffromIPPindex(obsS.station);
disp(['Calculate ' tag ' TEC'])
if isempty(refpos)
    try
        refpos = obsS.rcvpos(:)';
    catch
    end
end
% Geometry constants (shared); frequencies come from F
gpscons
Re = 6371.009*10^3;
h  = 350*10^3;
elev_mask = 30;
center_E = [0 0 0];
f1 = F.f1; f2 = F.f2; k = F.k;
lambda1 = F.lam1; lambda2 = F.lam2;

%% 1. Prepare matrix
NPRN = 64;
TEC.vertical         = nan(86400,NPRN);
TEC.slant            = nan(86400,NPRN);
TEC.withrcvbias      = nan(86400,NPRN);
TEC.withbias         = nan(86400,NPRN);
STECp                = nan(86400,NPRN);
STECl                = nan(86400,NPRN);
Times                = nan(86400,NPRN);
ROTI                 = nan(86400,NPRN);
DCB.sat              = nan(1,NPRN);
DCB.rcv              = nan;
prm.elevation        = nan(86400,NPRN);

%% 2. TEC calculation
Sat_obs = unique(obsS.index);
Sat_obs = Sat_obs(Sat_obs >= 1 & Sat_obs <= NPRN);
for i = 1:length(Sat_obs)
    disp(['PRN ' tag '#' num2str(Sat_obs(i)) ' ...'])
    try
    PRN = Sat_obs(i);
    Sat  = find(obsS.index == PRN);
    Time = round(obsS.epoch(Sat));
    if max(Time) > 86398
        Inx = Time+1 > 86399;
        Sat(Inx) = [];
        Time(Inx) = [];
    end
    if isempty(Sat), continue; end

    % 2.1 Read pseudorange / carrier phase (selected pair)
    C1 = obsS.data(Sat,ismember(obsS.type,PT.c1));
    P2 = obsS.data(Sat,ismember(obsS.type,PT.c2));
    L1 = lambda1*obsS.data(Sat,ismember(obsS.type,PT.l1));
    L2 = lambda2*obsS.data(Sat,ismember(obsS.type,PT.l2));

    % 2.2 elevation angle
    [satpos,~]  = satpos_xyz_sbias(Time,PRN,navS.eph,navS.index,C1);
    vector_s  = satpos-refpos;
    vector_r2 = refpos-center_E;
    vector_r  = repmat(vector_r2,length(vector_s),1);
    prm.elevation(Time+1,PRN) = 90-acosd(dot(vector_s,vector_r,2)./(vecnorm(vector_s')'...
                                .*vecnorm(vector_r')'));
    % 2.3 STEC
    STECp(Time+1,PRN) = k*(P2-C1);
    STECl(Time+1,PRN) = k*(L1-L2);
    Times(Time+1,PRN) = Time+1;
    catch
        disp(['PRN ' tag '#' num2str(Sat_obs(i)) '... error ...'])
        continue
    end
end

    % 2.4 elevation mask
mask                 = prm.elevation;
mask(mask<elev_mask) = NaN;
mask(~isnan(mask))   = 1;
TEC.STECp = mask.*STECp;
TEC.STECl = mask.*STECl;
prm.elevation = mask.*prm.elevation;
prm.Times = mask.*Times;

    % 2.5 Cycle slip correction (shared routine)
STECl_M_new = CycleSlipCorrection_v2(TEC.STECl,prm.elevation,prm.Times,Sat_obs');

    % 2.6 Carrier-to-code levelling per continuous arc
TEC.withbias = windowshiftTEC_SYS(STECl_M_new,TEC.STECp);
tec_min = nanmin(TEC.withbias(:));
if tec_min<0
    TEC.N_withbias    =  TEC.withbias + abs(tec_min);
else
    TEC.N_withbias    =  TEC.withbias;
end

    % 2.7 remove satellite bias (seconds -> TECU with pair k)
c  = 299792458;
Bias_sat_tec = satBias_s(:)'*(c*k);
DCB.sat      = Bias_sat_tec(1:NPRN);
disp('Remove satellite bias ....')
STEC_adj_nosatbias = TEC.N_withbias - ones(size(TEC.N_withbias))*diag(DCB.sat');
STEC_nooutline     = outlinecorr(STEC_adj_nosatbias);
TEC.withrcvbias    = STEC_nooutline;
for PNN = 1:NPRN
     count=0;
     for chk = 84900:86400
         if ~isnan(TEC.withrcvbias(chk,PNN))
             count = count+1;
         end
     end
     if count<1500
         TEC.withrcvbias(84900:86400,PNN)=nan;
     end
end

    % 2.8 receiver bias (Ma & Maruyama search, pair frequencies)
slant_factor       = sqrt(1-(Re*cosd(prm.elevation)/(Re+h)).^2);
disp('Remove receiver bias ....')
rcv_bias_ns = rcv_bias_sys(-30,30,0.1,TEC.withrcvbias,slant_factor,f1,f2);
A = 40.3;
DCB.rcv     = (rcv_bias_ns*10^-9*(c*(f1^2*f2^2/(A*(f1^2-f2^2)*10^16))));

STEC_completed = (TEC.withrcvbias - DCB.rcv);
VTEC_completed  = (TEC.withrcvbias - DCB.rcv).*slant_factor;
tec_min = nanmin(VTEC_completed(:));
TEC.slant    = STEC_completed - tec_min;
TEC.vertical = VTEC_completed - tec_min;

%% 3. ROTI (shared routine)
disp('Compute ROTI ....')
ROTI = roticalculation(STECl_M_new,Sat_obs');

%% 4. S4 amplitude scintillation index
% SNR selection: derive from the selected code pair (same signal as TEC/ROTI).
% PT.c1 = 'C1C' -> snrType = 'S1C' (C->S prefix substitution).
% Falls back to any S* on the same frequency if the exact match is absent.
disp('Compute S4 amplitude scintillation ....')
snrType_sys = ['S' PT.c1(2:end)];   % e.g. C1C->S1C, B1C->S1C, E1C->S1C
SNR_col_idx_sys = find(ismember(obsS.type, snrType_sys), 1);
if isempty(SNR_col_idx_sys)
    % Fallback: any S1* or S* column present for this constellation
    s_cands = obsS.type(strncmp(obsS.type, 'S1', 2));
    if isempty(s_cands)
        s_cands = obsS.type(strncmp(obsS.type, 'S',  1));
    end
    if ~isempty(s_cands)
        snrType_sys    = s_cands{1};
        SNR_col_idx_sys = find(ismember(obsS.type, snrType_sys), 1);
        disp([tag ' S4: SNR fallback to ' snrType_sys])
    end
end

S4_SYS = nan(86400, NPRN);
if ~isempty(SNR_col_idx_sys)
    SNR_mat_sys = nan(86400, NPRN);
    for i2 = 1:length(Sat_obs)
        try
            PRN2  = Sat_obs(i2);
            Sat2  = find(obsS.index == PRN2);
            Time2 = round(obsS.epoch(Sat2));
            if max(Time2) > 86398
                Inx2  = Time2+1 > 86399;
                Sat2(Inx2)  = [];
                Time2(Inx2) = [];
            end
            if isempty(Sat2), continue; end
            snr_raw = obsS.data(Sat2, SNR_col_idx_sys);
            % Apply elevation mask
            elev_ok = ~isnan(prm.elevation(Time2+1, PRN2));
            snr_raw(~elev_ok) = NaN;
            SNR_mat_sys(Time2+1, PRN2) = snr_raw;
        catch
        end
    end
    S4_SYS = S4calculation(SNR_mat_sys, Sat_obs');
    disp([tag ' S4: ' snrType_sys ' (' num2str(sum(isfinite(S4_SYS(:)))) ' valid epochs)'])
else
    disp([tag ' S4: no SNR column found — S4 matrix all-NaN'])
end

%% Save file
year  = num2str(obsS.date(1));
month = num2str(obsS.date(2),'%.2d');
date  = num2str(obsS.date(3),'%.2d');
name1 = ['TEC_' tag '_' year '_' month '_' date];
name2 = ['DCB_' tag '_' year '_' month '_' date];
name3 = ['ROTI_' tag '_' year '_' month '_' date];
name4 = ['prm_' tag '_' year '_' month '_' date];
name5 = ['S4_' tag '_' year '_' month '_' date];
eval([name1 '= TEC;'])
eval([name2 '= DCB;'])
eval([name3 '= ROTI;'])
eval([name4 '= prm;'])
eval([name5 '= S4_SYS;'])
filename = [S_path 'TEC_' tag '_' year '_' month '_' date];
save(filename,name1,name2,name3,name4,name5,'refpos')
disp(['Complete ' tag ' TEC'])
end

function STEC_adj = windowshiftTEC_SYS(STECl,STECp)
% windowshiftTEC generalized to any PRN count (same algorithm as original)
nPRN = size(STECl, 2);
STEC_adj = nan(size(STECl));
for LP = 1:nPRN
    [val,~] = find(~isnan(STECl(:,LP)));
    if isempty(val), continue; end
    zz = find(diff(val)>=10800); %3hr*3600
    if isempty(zz)
        STEC_adj(:,LP) = STECl(:,LP)+ ones(size(STECl(:,LP)))...
          *diag(nanmean(STECp(:,LP)-STECl(:,LP)));
        continue
    end
    wind = [];
    for wd = 1:length(zz)
        if val(zz(wd)) <=82800
            wind(wd) = val(zz(wd))+3600; %#ok<AGROW>
        end
    end
    if isempty(wind)
        STEC_adj(:,LP) = STECl(:,LP)+ ones(size(STECl(:,LP)))...
          *diag(nanmean(STECp(:,LP)-STECl(:,LP)));
        continue
    end
    wind(wd+1) = length(STECl);
    for SW = 1:length(wind)
        if SW == 1
            STEC_adj(1:wind(SW),LP)=STECl(1:wind(SW),LP) ...
            + ones(size(STECl(1:wind(SW),LP)))...
            *diag(nanmean(STECp(1:wind(SW),LP)...
            -STECl(1:wind(SW),LP)));
        else
            STEC_adj(wind(SW-1):wind(SW),LP)=STECl(wind(SW-1):wind(SW),LP) ...
            + ones(size(STECl(wind(SW-1):wind(SW),LP)))...
            *diag(nanmean(STECp(wind(SW-1):wind(SW),LP)...
            -STECl(wind(SW-1):wind(SW),LP)));
        end
    end
end
end
