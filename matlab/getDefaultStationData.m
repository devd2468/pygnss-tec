function stationData = getDefaultStationData()
% Default station positions (fallback)

    stations = {
        'RUTI', [-1293985, 6025319, 1638960];
        'LCK4', [-1294882, 6024986, 1639293];
        'GDKG', [-1201408, 6066337, 1555734];
        'BHPL', [-1140009, 6061593, 1618814];
        'IISC', [-1158290, 6087925, 1503753];
        'DRDN', [-1176791, 6096212, 1455314];
        'IITK', [-1326863, 5923373, 1951940];
        'SHLG', [-1078589, 6046539, 1713801];
        'PBR4', [-1146424, 6089932, 1504580];
        'JDPR', [-1294882, 6024986, 1639293];
        'SRTN', [-1021155, 6214494, 1005605];
        'KMIT', [-1158290, 6087925, 1503753];
    };
    stationData = repmat(struct('name', '', 'xyz', []), length(stations), 1);
    for i = 1:length(stations)
        stationData(i).name = stations{i, 1};
        stationData(i).xyz = stations{i, 2};
    end
end
