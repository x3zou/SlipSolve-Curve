function startup(repoRoot)
%STARTUP Add SlipSolve source code to the MATLAB path.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
end

addpath(fullfile(repoRoot, "src"));

end

