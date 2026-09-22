function lla = ecef2lla(XYZ)
% Convert ECEF coordinates to geodetic latitude, longitude, altitude

    X = XYZ(1); Y = XYZ(2); Z = XYZ(3);
    a = 6378137.0;
    f = 1/298.257223563;
    e2 = 2*f - f^2;
    lon = atan2(Y, X);
    p = sqrt(X^2 + Y^2);
    lat = atan2(Z, p*(1-e2));
    for k = 1:10
        N = a / sqrt(1 - e2*sin(lat)^2);
        alt = p/cos(lat) - N;
        lat_new = atan2(Z + e2*N*sin(lat), p);
        if abs(lat_new - lat) < 1e-12
            lat = lat_new; break;
        end
        lat = lat_new;
    end
    lla = [lat*180/pi, lon*180/pi, alt];
end
