function S = aloadTEC(dayRoot, station4)
% ALOADTEC  Load one station-day MultiGNSS result struct.
% Returns [] if missing. Fields used: TEC, ROTI, DCB, prm, station.
    S = [];
    d = dir(fullfile(dayRoot, 'Results', sprintf('MultiGNSS_%s_*.mat', station4)));
    if isempty(d)
        return;
    end
    try
        S = load(fullfile(d(1).folder, d(1).name));
    catch
        S = [];
    end
end
