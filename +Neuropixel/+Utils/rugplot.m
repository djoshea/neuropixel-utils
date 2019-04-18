function h = rugplot(ticks, varargin)

    if nargin == 0
        pos = randn(500, 4);
        cla;
        xlim([-4 4]);
        ylim([-4 4]);
        hold on;
        Neuropixel.Utils.rugplot(pos(:, 1), 'side', 'top', 'Color', 'r');
        Neuropixel.Utils.rugplot(pos(:, 2), 'side', 'right', 'Color', 'g');
        Neuropixel.Utils.rugplot(pos(:, 3), 'side', 'bottom', 'Color', 'c');
        Neuropixel.Utils.rugplot(pos(:, 4), 'side', 'left', 'Color', 'm');
        return;
    end

    p = inputParser();
    p.addParameter('Color', [0 0 1], @isvector);
    p.addParameter('side', 'top', @ischar);
    p.addParameter('length', 0.01, @isscalar); % fraction of axis extents
    p.addParameter('offset', 0, @isscalar); % fraction of axis extents, + mean away from center of axis
    p.addParameter('expand_limits', false, @islogical);
    p.parse(varargin{:});
    
    side = p.Results.side;
    ylfull = ylim;
    xlfull = xlim;  
    expand = p.Results.expand_limits;

    switch side
        case {'top', 'bottom'}
            full = ylfull;
        case {'left', 'right'}
            full = xlfull;
        otherwise
            error('Unknown side');
    end
    len = p.Results.length * (full(2)-full(1));
    
    offset = p.Results.offset;
       
    switch side
        case 'top'
            if expand
                ylfull(2) = ylfull(2) + len;
                ylim(ylfull);
            end
            x = [ticks'; ticks'];
            y = repmat([ylfull(2); ylfull(2)-len], 1, numel(ticks)) + offset;
        case 'bottom'
            if expand
                ylfull(1) = ylfull(1) - len;
                ylim(ylfull);
            end
            x = [ticks'; ticks'];
            y = repmat([ylfull(1); ylfull(1)+len], 1, numel(ticks)) - offset;
        case 'right'
            if expand
                xlfull(2) = xlfull(2) + len;
                xlim(xlfull);
            end
            y = [ticks'; ticks'];
            x = repmat([xlfull(2); xlfull(2)-len], 1, numel(ticks)) + offset;
        case 'left'
            if expand
                xlfull(1) = ylfull(1) - len;
                xlim(xlfull);
            end
            y = [ticks'; ticks'];
            x = repmat([xlfull(1); xlfull(1)+len], 1, numel(ticks)) - offset;
    end
    hold on;
    h = line(x, y, 'Color', p.Results.Color, 'LineWidth', 0.5);
end