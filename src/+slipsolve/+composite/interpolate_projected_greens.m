function green = interpolate_projected_greens(componentFiles, xVector, yVector, ...
    x, y, projection, nColumns, settings)
%INTERPOLATE_PROJECTED_GREENS Sample full XYZ Greens and project observations.

xVector = double(xVector(:).');
yVector = double(yVector(:).');
x = double(x(:));
y = double(y(:));
projection = double(projection);
if numel(x) ~= numel(y) || size(projection, 1) ~= numel(x) || ...
        size(projection, 2) ~= 3
    error("SlipSolve:InvalidInterpolationInput", ...
        "x, y, and the N-by-3 projection matrix must have matching rows.");
end
if any(~isfinite([x; y]), "all") || any(~isfinite(projection), "all")
    error("SlipSolve:InvalidInterpolationInput", ...
        "Coordinates and projection coefficients must be finite.");
end
if any(x < min(xVector) | x > max(xVector) | ...
        y < min(yVector) | y > max(yVector))
    error("SlipSolve:InterpolationOutsideGrid", ...
        "One or more observations lie outside the full composite grid.");
end

method = string(get_field(settings, "method", "linear"));
if ~strcmpi(method, "linear")
    error("SlipSolve:InvalidInterpolationConfig", ...
        "Full-grid interpolation currently supports method='linear'.");
end
blockSize = double(get_field(settings, "columnBlockSize", 32));
if ~isscalar(blockSize) || blockSize < 1 || mod(blockSize, 1) ~= 0
    error("SlipSolve:InvalidInterpolationConfig", ...
        "columnBlockSize must be a positive integer.");
end

components = [ ...
    struct("file", string(componentFiles.east), ...
           "variable", "G_e", "coefficient", 1); ...
    struct("file", string(componentFiles.north), ...
           "variable", "G_n", "coefficient", 2); ...
    struct("file", string(componentFiles.vertical), ...
           "variable", "G_u", "coefficient", 3)];
expectedRows = numel(xVector)*numel(yVector);
interpolationMatrix = bilinear_interpolation_matrix(xVector, yVector, x, y);
green = zeros(numel(x), nColumns);

for c = 1:numel(components)
    coefficient = projection(:, components(c).coefficient);
    if ~any(coefficient)
        continue
    end
    validate_component(components(c).file, components(c).variable, ...
        expectedRows, nColumns);
    fprintf("Interpolating %s from full grid to %d observations.\n", ...
        components(c).variable, numel(x));
    loaded = load(components(c).file, components(c).variable);
    sourceValues = loaded.(components(c).variable);
    clear loaded
    for first = 1:blockSize:nColumns
        last = min(first+blockSize-1, nColumns);
        columns = first:last;
        values = interpolationMatrix*double(sourceValues(:, columns));
        green(:, columns) = green(:, columns) + coefficient.*values;
        fprintf("  %s columns %d-%d / %d\n", ...
            components(c).variable, first, last, nColumns);
    end
    clear sourceValues values
end
if any(~isfinite(green), "all")
    error("SlipSolve:InterpolationFailed", ...
        "Projected sampled composite Green contains NaN or Inf values.");
end
end

function matrix = bilinear_interpolation_matrix(xVector, yVector, x, y)
nx = numel(xVector);
ny = numel(yVector);
if nx < 2 || ny < 2 || any(diff(xVector) <= 0) || any(diff(yVector) <= 0)
    error("SlipSolve:InvalidCompositeGrid", ...
        "Composite x/y vectors must be strictly increasing with at least two values.");
end
xEdges = [-Inf, xVector(2:end-1), Inf];
yEdges = [-Inf, yVector(2:end-1), Inf];
ix = discretize(x, xEdges);
iy = discretize(y, yEdges);
ix = min(max(ix, 1), nx-1);
iy = min(max(iy, 1), ny-1);
tx = (x-xVector(ix).') ./ (xVector(ix+1).'-xVector(ix).');
ty = (y-yVector(iy).') ./ (yVector(iy+1).'-yVector(iy).');
row = (1:numel(x)).';
i00 = sub2ind([ny nx], iy, ix);
i10 = sub2ind([ny nx], iy, ix+1);
i01 = sub2ind([ny nx], iy+1, ix);
i11 = sub2ind([ny nx], iy+1, ix+1);
matrix = sparse([row; row; row; row], [i00; i10; i01; i11], ...
    [(1-tx).*(1-ty); tx.*(1-ty); (1-tx).*ty; tx.*ty], ...
    numel(x), nx*ny);
end

function validate_component(filePath, variableName, expectedRows, expectedColumns)
if exist(filePath, "file") ~= 2
    error("SlipSolve:MissingInput", "Composite component is missing: %s", filePath);
end
info = whos("-file", filePath, variableName);
if isempty(info) || ~isequal(info.size, [expectedRows expectedColumns])
    error("SlipSolve:GreenSizeMismatch", ...
        "%s must have size [%d %d] in %s.", ...
        variableName, expectedRows, expectedColumns, filePath);
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
