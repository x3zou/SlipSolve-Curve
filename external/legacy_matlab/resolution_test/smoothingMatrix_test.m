function [H, h1] = smoothingMatrix_test(n1,varargin)
% smoothingMatrix constructs a Tikhonov smoothing matrix H that
% penalizes the difference in slip between neighboring fault patches.
% Xiaoyu Zou, 04/06/2025
% Input:
%   n1 - A (np, X) contact matrix. Each row corresponds to a slip patch,
%        and the three columns give the indices of neighboring patches.
%        A NaN value indicates that there is no neighboring patch in that spot.
%
% Output:
%   H  - A (2*numConstraints x (2*np)) sparse smoothing matrix.
%        The first np columns correspond to strike-slip, and the next np
%        columns correspond to dip-slip.
%   h1 - size of H
%
% The smoothing constraints are of the form:
%   For strike-slip:  m(i) - m(neighbor) = 0
%   For dip-slip:     m(np+i) - m(np+neighbor) = 0

% Alternatively, it can constrin the difference between the central element
% and the averaged slip on the neibors:
%   For strike-slip:  m(i) - m(neighbor)/N = 0
%   For dip-slip:     m(np+i) - m(np+neighbor)/N = 0
r1=1;% relative smoothness for strike-slip
r2=1;% relative smoothness for dip-slip
horizontal_smooth_layer=[];
vertical_smooth_layer=[];
if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'ss_ratio'
                    r1=varargin{CC*2};
                case 'ds_ratio'
                    r2=varargin{CC*2};
                case 'layer'
                    layer=varargin{CC*2};
                case 'horizontal_smooth_layer'
                    horizontal_smooth_layer=varargin{CC*2};
                case 'horizontal_smooth_ratio'
                    horizontal_smooth_ratio=varargin{CC*2};
                case 'vertical_smooth_layer'
                    vertical_smooth_layer=varargin{CC*2};
                case 'vertical_smooth_ratio'
                    vertical_smooth_ratio=varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

    % Number of slip patches
    np = size(n1, 1);
    % Total number of model parameters (strike-slip and dip-slip)
    N = 2 * np;
    
    % Count the total number of valid neighbor constraints
    numConstraints = sum(~isnan(n1(:)));
    % Each neighbor yields two constraints (one for strike, one for dip)
    totalRows = 2 * numConstraints;
    
    % Preallocate arrays for building the sparse matrix L.
    % Each constraint contributes two non-zero entries (+1 and -1).
    numEntries = 2 * totalRows;
    rows = zeros(numEntries, 1);
    cols = zeros(numEntries, 1);
    vals = zeros(numEntries, 1);
    
    rowCounter = 0;     % Keeps track of the current row in L
    entryCounter = 0;   % Keeps track of the current index in rows/cols/vals arrays
    
    % Loop over each slip patch and its potential neighbors
    for i = 1:np
        for j = 1:size(n1,2)
            neighbor = n1(i,j);
            valid_neighbor=sum(~isnan(n1(i,:))); % added for test
            if ~isnan(neighbor)
                % For strike-slip: Add a constraint row (m(i) - m(neighbor) = 0)
                rowCounter = rowCounter + 1;
                entryCounter = entryCounter + 1;
                rows(entryCounter) = rowCounter;
                cols(entryCounter) = i;  % strike-slip for patch i
                vals(entryCounter) = r1;
                if ~isempty(horizontal_smooth_layer)
                    if any(horizontal_smooth_layer==layer(i)) && layer(i) == layer(neighbor)
                        vals(entryCounter) = r1*horizontal_smooth_ratio;
                    end
                end

                if ~isempty(vertical_smooth_layer)
                    if any(vertical_smooth_layer==layer(i)) && layer(i) ~=layer(neighbor)
                        vals(entryCounter) = r1*vertical_smooth_ratio;
                    end
                end
                
                entryCounter = entryCounter + 1;
                rows(entryCounter) = rowCounter;
                cols(entryCounter) = neighbor;  % strike-slip for neighbor patch
                vals(entryCounter) = -1*r1; %default
                if ~isempty(horizontal_smooth_layer)
                    if any(horizontal_smooth_layer==layer(i)) && layer(i) == layer(neighbor)
                        vals(entryCounter) = -1*r1*horizontal_smooth_ratio;
                    end
                end

                if ~isempty(vertical_smooth_layer)
                    if any(vertical_smooth_layer==layer(i)) && layer(i) ~=layer(neighbor)
                        vals(entryCounter) = -1*r1*vertical_smooth_ratio;
                    end
                end
                %vals(entryCounter) = -1/valid_neighbor; % averaged
                
                % For dip-slip: Add a constraint row (m(np+i) - m(np+neighbor) = 0)
                rowCounter = rowCounter + 1;
                entryCounter = entryCounter + 1;
                rows(entryCounter) = rowCounter;
                cols(entryCounter) = np + i;  % dip-slip for patch i
                vals(entryCounter) = r2;
                if ~isempty(horizontal_smooth_layer)
                    if any(horizontal_smooth_layer==layer(i)) && layer(i) == layer(neighbor)
                        vals(entryCounter) = r2*horizontal_smooth_ratio;
                    end
                end

                if ~isempty(vertical_smooth_layer)
                    if any(vertical_smooth_layer==layer(i)) && layer(i) ~=layer(neighbor)
                        vals(entryCounter) = r2*vertical_smooth_ratio;
                    end
                end
                
                entryCounter = entryCounter + 1;
                rows(entryCounter) = rowCounter;
                cols(entryCounter) = np + neighbor;  % dip-slip for neighbor patch
                vals(entryCounter) = -1*r2; %default
                if ~isempty(horizontal_smooth_layer)
                    if any(horizontal_smooth_layer==layer(i)) && layer(i) == layer(neighbor)
                        vals(entryCounter) = -1*r2*horizontal_smooth_ratio;
                    end
                end

                if ~isempty(vertical_smooth_layer)
                    if any(vertical_smooth_layer==layer(i)) && layer(i) ~=layer(neighbor)
                        vals(entryCounter) = -1*r2*vertical_smooth_ratio;
                    end
                end
                %vals(entryCounter) = -1/valid_neighbor; %averaged
            end
        end
    end
    
    % Build the sparse smoothing matrix L
    H = sparse(rows, cols, vals, totalRows, N);
    h1 = size(H,1);
end