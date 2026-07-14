function report = verify_example_data(repoRoot)
%VERIFY_EXAMPLE_DATA Check the downloaded Google Drive example-data bundle.

if nargin < 1 || strlength(string(repoRoot)) == 0
    exampleDirectory = fileparts(mfilename("fullpath"));
    repoRoot = fileparts(fileparts(exampleDirectory));
end
repoRoot = string(repoRoot);
manifestFile = fullfile(fileparts(mfilename("fullpath")), ...
    "example_data_manifest.csv");
manifest = readtable(manifestFile, "TextType", "string");

status = strings(height(manifest), 1);
actualBytes = nan(height(manifest), 1);
for k = 1:height(manifest)
    relativePath = strrep(manifest.relative_path(k), "/", filesep);
    target = fullfile(repoRoot, relativePath);
    info = dir(target);
    if isempty(info)
        status(k) = "missing";
        continue
    end
    actualBytes(k) = info(1).bytes;
    if actualBytes(k) == manifest.bytes(k)
        status(k) = "ok";
    else
        status(k) = "wrong_size";
    end
end

report = manifest;
report.actual_bytes = actualBytes;
report.status = status;
bad = status ~= "ok";
if any(bad)
    details = join(report.relative_path(bad)+" ["+status(bad)+"]", newline);
    error("SlipSolveCurve:IncompleteExampleData", ...
        "Google Drive example data are missing or incomplete:\n%s\n\n"+ ...
        "Download SlipSolve-curve-example-data and copy its data "+ ...
        "and greens folders into the repository root.", details);
end

totalGiB = sum(double(manifest.bytes)) / 1024^3;
fprintf("Verified %d Google Drive example-data files (%.2f GiB).\n", ...
    height(manifest), totalGiB);
end
