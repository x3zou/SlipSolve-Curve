function cfg = load(configPath)
%LOAD Load a SlipSolve configuration function.

if nargin < 1 || strlength(string(configPath)) == 0
    configPath = fullfile(pwd, "config", "example_project.m");
end

configPath = char(configPath);
if exist(configPath, "file") ~= 2
    error("SlipSolve:ConfigNotFound", "Config file not found: %s", configPath);
end

[configFolder, configName] = fileparts(configPath);
previousFolder = pwd;
cleanup = onCleanup(@() cd(previousFolder));
cd(configFolder);

configFunction = str2func(configName);
cfg = configFunction();

if ~isstruct(cfg)
    error("SlipSolve:InvalidConfig", "Config function must return a struct.");
end

end

