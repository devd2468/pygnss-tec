function r = repoRoot()
% REPOROOT  Directory containing this repository's code.
% Portable alternative to pwd-based lookups: works no matter where MATLAB
% was launched from, on any OS (uses filesep-safe fullfile at call sites).
r = fileparts(mfilename('fullpath'));
if isempty(r)
    r = pwd;
end
end
