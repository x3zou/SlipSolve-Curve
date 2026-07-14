function products = convert_legacy_getedgrn(inputFile, outputDirectory, overwrite)
%CONVERT_LEGACY_GETEDGRN Convert EDGRN text output using getedgrn.m logic.
%
% This is a path-safe copy of the legacy getedgrn.m algorithm. It reads the
% three output-file names from the fifth non-comment EDGRN input line, reads
% each database with textread, and saves the same variables as legacy files.

if nargin < 3
    overwrite = false;
end

inputFile = string(inputFile);
outputDirectory = string(outputDirectory);
if exist(inputFile, "file") ~= 2
    error("SlipSolve:MissingInput", "EDGRN input file does not exist: %s", inputFile);
end
if exist(outputDirectory, "dir") ~= 7
    mkdir(outputDirectory);
end

[databaseDirectory, databaseNames] = read_edgrn_output_names(inputFile);
[~, inputBase, inputExtension] = fileparts(inputFile);
inputLabel = inputBase + inputExtension;
types = ["ss", "ds", "cl"];
products = repmat(struct( ...
    "type", "", "textFile", "", "matFile", "", "nr", 0, "nz", 0, ...
    "r1", 0, "r2", 0, "z1", 0, "z2", 0, "receiverDepth", 0), 3, 1);

for k = 1:numel(types)
    textFile = fullfile(databaseDirectory, databaseNames(k));
    if exist(textFile, "file") ~= 2
        error("SlipSolve:MissingInput", ...
            "EDGRN output declared by %s is missing: %s", inputFile, textFile);
    end

    matFile = fullfile(outputDirectory, "edgrn_" + types(k) + "_" + inputLabel + ".mat");
    if exist(matFile, "file") == 2 && ~overwrite
        error("SlipSolve:ExistingOutput", ...
            "Converted EDGRN file already exists: %s. Set overwriteExistingMatFiles=true to replace it.", matFile);
    end

    database = read_fundamental_database(textFile);
    nr = database.nr;
    nz = database.nz;
    r1 = database.r1;
    r2 = database.r2;
    z1 = database.z1;
    z2 = database.z2;
    zrec0 = database.zrec0;
    lambda = database.lambda;
    mu = database.mu;
    uz = database.uz;
    ur = database.ur;
    ut = database.ut;
    ezz = database.ezz;
    err = database.err;
    ett = database.ett;
    ezr = database.ezr;
    ert = database.ert;
    etz = database.etz;
    duz_dr = database.duz_dr;

    % Variable names and source reader intentionally match legacy getedgrn.m.
    save(matFile, "nr", "nz", "r1", "r2", "z1", "z2", "zrec0", "lambda", "mu", ...
        "uz", "ur", "ut", "ezz", "err", "ett", "ezr", "ert", "etz", "duz_dr");

    products(k).type = types(k);
    products(k).textFile = textFile;
    products(k).matFile = matFile;
    products(k).nr = nr;
    products(k).nz = nz;
    products(k).r1 = r1;
    products(k).r2 = r2;
    products(k).z1 = z1;
    products(k).z2 = z2;
    products(k).receiverDepth = zrec0;
end
end

function [databaseDirectory, names] = read_edgrn_output_names(inputFile)
fid = fopen(inputFile, "r");
if fid < 0
    error("SlipSolve:FileOpen", "Could not open EDGRN input file: %s", inputFile);
end
cleanup = onCleanup(@() fclose(fid));

nonCommentLine = 0;
outputLine = "";
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~contains(string(line), "#")
        nonCommentLine = nonCommentLine + 1;
        if nonCommentLine == 5
            outputLine = string(line);
            break
        end
    end
end

quoted = regexp(outputLine, "'([^']*)'", "tokens");
if numel(quoted) < 4
    error("SlipSolve:InvalidEdgrnInput", ...
        "Could not read EDGRN output directory and SS/DS/CL names from the fifth non-comment line of %s.", inputFile);
end

rawDirectory = string(quoted{1}{1});
names = [string(quoted{2}{1}), string(quoted{3}{1}), string(quoted{4}{1})];
[inputDirectory, ~] = fileparts(inputFile);
if startsWith(rawDirectory, "/") || ~isempty(regexp(rawDirectory, "^[A-Za-z]:[\\\\/]", "once"))
    databaseDirectory = rawDirectory;
else
    databaseDirectory = fullfile(inputDirectory, rawDirectory);
end
end

function database = read_fundamental_database(textFile)
fid = fopen(textFile, "r");
if fid < 0
    error("SlipSolve:FileOpen", "Could not open EDGRN output file: %s", textFile);
end
cleanup = onCleanup(@() fclose(fid));

headerLines = 0;
metadataLine = "";
while ~feof(fid)
    line = fgetl(fid);
    headerLines = headerLines + 1;
    if ischar(line) && ~contains(string(line), "#")
        metadataLine = string(line);
        break
    end
end

metadata = sscanf(strrep(strrep(metadataLine, "D", "E"), "d", "e"), "%f");
if numel(metadata) < 9
    error("SlipSolve:InvalidEdgrnOutput", "Invalid EDGRN database header in %s.", textFile);
end

% textread is deliberately retained to reproduce the legacy getedgrn.m reader.
[uz, ur, ut, ezz, err, ett, ezr, ert, etz, duz_dr] = textread( ...
    char(textFile), "%f %f %f %f %f %f %f %f %f %f\n ", "headerlines", headerLines); %#ok<DTXTRD>

expectedRows = metadata(1) * metadata(4);
if numel(uz) ~= expectedRows
    error("SlipSolve:InvalidEdgrnOutput", ...
        "EDGRN database %s has %d rows; expected nr*nz = %d.", ...
        textFile, numel(uz), expectedRows);
end

database = struct();
database.nr = metadata(1);
database.r1 = metadata(2);
database.r2 = metadata(3);
database.nz = metadata(4);
database.z1 = metadata(5);
database.z2 = metadata(6);
database.zrec0 = metadata(7);
database.lambda = metadata(8);
database.mu = metadata(9);
database.uz = uz;
database.ur = ur;
database.ut = ut;
database.ezz = ezz;
database.err = err;
database.ett = ett;
database.ezr = ezr;
database.ert = ert;
database.etz = etz;
database.duz_dr = duz_dr;
end
