function tracks = read_sample_list(sampleListFile)
%READ_SAMPLE_LIST Read a legacy InSAR sample-list text file.
%
% The current initial sampler expects fourteen whitespace-delimited columns:
% track directory, grid file, data type, nominal npt, lon/lat bounds,
% min/max cell size in km, RMS threshold, output MAT name, and x/y padding.

if nargin < 1 || strlength(string(sampleListFile)) == 0
    error("SlipSolve:MissingInput", "A sample-list file path is required.");
end

sampleListFile = char(sampleListFile);
if exist(sampleListFile, "file") ~= 2
    error("SlipSolve:MissingInput", "Sample-list file not found: %s", sampleListFile);
end

lines = readlines(sampleListFile);
lines = strip(lines);
lines(lines == "" | startsWith(lines, "#")) = [];

tracks = repmat(empty_track(), numel(lines), 1);
baseDir = string(fileparts(sampleListFile));

for k = 1:numel(lines)
    parts = split(lines(k));
    parts(parts == "") = [];

    if numel(parts) < 12
        error("SlipSolve:InvalidSampleList", ...
            "Sample-list row %d has %d columns; expected at least 12.", k, numel(parts));
    end

    track = empty_track();
    track.row = k;
    track.sourceLine = lines(k);
    track.directory = string(parts(1));
    track.directoryFullPath = normalize_path(baseDir, track.directory);
    track.gridFile = string(parts(2));
    track.gridFullPath = fullfile(track.directoryFullPath, track.gridFile);
    track.type = string(parts(3));
    track.nominalPointCount = str2double(parts(4));
    track.boundsLonLat = str2double(parts(5:8)).';
    track.minCellSizeKm = str2double(parts(9));
    track.maxCellSizeKm = str2double(parts(10));
    track.rmsThreshold = str2double(parts(11));
    track.outputFile = string(parts(12));
    track.outputFullPath = fullfile(track.directoryFullPath, track.outputFile);

    if numel(parts) >= 14
        track.padX = str2double(parts(13));
        track.padY = str2double(parts(14));
    end

    validate_track(track, k);
    tracks(k) = track;
end

end

function track = empty_track()
track = struct();
track.row = NaN;
track.sourceLine = "";
track.directory = "";
track.directoryFullPath = "";
track.gridFile = "";
track.gridFullPath = "";
track.type = "";
track.nominalPointCount = NaN;
track.boundsLonLat = nan(1, 4);
track.minCellSizeKm = NaN;
track.maxCellSizeKm = NaN;
track.rmsThreshold = NaN;
track.outputFile = "";
track.outputFullPath = "";
track.padX = 0;
track.padY = 0;
end

function pathOut = normalize_path(baseDir, pathIn)
pathIn = string(pathIn);
if startsWith(pathIn, filesep) || startsWith(pathIn, "~")
    pathOut = pathIn;
else
    pathOut = fullfile(baseDir, pathIn);
end
end

function validate_track(track, row)
numericValues = [track.nominalPointCount, track.boundsLonLat, ...
    track.minCellSizeKm, track.maxCellSizeKm, track.rmsThreshold, ...
    track.padX, track.padY];

if any(isnan(numericValues))
    error("SlipSolve:InvalidSampleList", "Sample-list row %d contains a nonnumeric required value.", row);
end

if track.minCellSizeKm <= 0 || track.maxCellSizeKm <= 0
    error("SlipSolve:InvalidSampleList", "Sample-list row %d has nonpositive cell-size bounds.", row);
end

if track.minCellSizeKm > track.maxCellSizeKm
    error("SlipSolve:InvalidSampleList", "Sample-list row %d has min cell size greater than max cell size.", row);
end
end

