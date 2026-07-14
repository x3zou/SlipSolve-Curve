function save_stage(path, variableName, value)
%SAVE_STAGE Save one named stage output in a MAT file.

folder = fileparts(path);
if exist(folder, "dir") ~= 7
    mkdir(folder);
end

payload = struct();
payload.(variableName) = value;
save(path, "-struct", "payload", "-v7.3");

end

