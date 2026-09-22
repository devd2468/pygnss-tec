function [dcb_s, nHit] = osbPairDCB(osb, sysChar, prnList, code1, code2)
% OSBPAIRDCB  Satellite DCB (code1-code2 pair) in seconds from OSB struct.
% Uses DCB(code1-code2) = OSB_code2 - OSB_code1  [seconds].
% Output dcb_s is 64x1 indexed by PRN (0 where unavailable); nHit = count.
    dcb_s = zeros(64, 1);
    nHit = 0;
    if isempty(osb)
        return;
    end
    for i = 1:numel(prnList)
        p = prnList(i);
        if p < 1 || p > 64, continue; end
        f = sprintf('%s%02d', sysChar, p);
        if isfield(osb, f)
            s = osb.(f);
            if isfield(s, code1) && isfield(s, code2)
                dcb_s(p) = (s.(code2) - s.(code1)) * 1e-9;
                nHit = nHit + 1;
            end
        end
    end
end
