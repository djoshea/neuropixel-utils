function map = convertColorsToMatrix(colors)

    hasAlpha = false;

    if ischar(colors)
        map = convStr(colors);
    elseif iscell(colors)
        N = numel(colors);
        map = nan(N, 4);

        for i = 1:N
            if ischar(colors{i})
                map(i, 1:3) = convStr(colors{i});
            elseif isvector(colors{i})
                n = numel(colors{i});
                if n == 3
                    map(i, 1:3) = colors{i};
                elseif n == 4
                    map(i, :) = colors{i};
                    hasAlpha = true;
                else
                    error('Numeric cell contents must be 3 or 4 vectors');
                end
            else
                error('Cell contents must be char or vectors');
            end
        end

        if ~hasAlpha
            map = map(:, 1:3);
        end

    elseif ismatrix(colors)
        assert(size(colors, 2) == 3 || size(colors, 2) == 4, 'Color matrix must have 3 or 4 columns'); 
        map = colors;
    else
        error('Colors must be cell array or matrix');
    end
end

function c = convStr(s)
    switch s
        case 'none'
            c = 'none';
        case 'k'
            c = [0 0 0];
        case 'b'
            c = [0 0 1];
        case 'g'
            c = [0 1 0];
        case 'c'
            c = [0 1 1];
        case 'r'
            c = [1 0 0];
        case 'm'
            c = [1 0 1];
        case 'y'
            c = [1 1 0];
        case 'w'
            c = [1 1 1];
        otherwise
            error('Unknown color %s', s);
    end
end