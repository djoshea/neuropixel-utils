function hs = errorshadeInterval(xmat, lomat, himat, color, varargin)
% hs = errorshadeInterval(x, lo, hi, color, varargin)

    p = inputParser();
    p.addParameter('errorColor', [], @(x) true);
    p.addParameter('shadeArgs', {}, @iscell);
    p.addParameter('axh', [], @(x) true);
    p.addParameter('alpha', 1, @isscalar);
    p.addParameter('z', 0, @isscalar); % used for visual stacking on 2-d plots
    p.addParameter('Clipping', true, @islogical);
    p.parse(varargin{:}); 
    
    color = Neuropixel.Utils.convertColorsToMatrix(color);
    
    z = p.Results.z;

    if isempty(p.Results.axh)
        axh = newplot;
    else
        axh = p.Results.axh;
    end

    cols = size(xmat, 2);
    hs = gobjects(cols, 1);

    for iC = 1:cols
        x = xmat(:, iC)';
        y1 = lomat(:, iC)';
        y2 = himat(:, iC)';
        
        if all(isnan(y1) | isnan(y2))   
            return;
        end
    
        
    % 
    %     % plot the shaded area
    %     x = makerow(x);
    %     y1 = makerow(y1);
    %     y2 = makerow(y2);
    
        % desaturate the color for shading if not translucent
        if isempty(p.Results.errorColor)
            shadeColor = Neuropixel.Utils.brighten(color, 0.5);
        else
            shadeColor = p.Results.errorColor;
        end
        
        % need to split the vecs
        nanMask = isnan(y1) | isnan(y2);
        offset = 1;
        while(offset < numel(x))
            % find next non-nan sample
            newOffset = find(~nanMask(offset:end), 1, 'first');
            if isempty(newOffset)
                break;
            end
            
            offset = offset+newOffset-1;
            nextNaN = find(nanMask(offset:end), 1, 'first');
            if isempty(nextNaN)
                regionEnd = numel(x);
            else
                regionEnd = nextNaN+offset - 2;
            end
            
            regionStart = offset;
            mask = regionStart:regionEnd;
            
            hs(iC) = shadeSimple(axh, x(mask), y1(mask), y2(mask), z, 'FaceColor', shadeColor, ...
                'alpha', p.Results.alpha, 'Clipping', p.Results.Clipping, p.Results.shadeArgs{:});

            hold(axh, 'on');
            
            offset = regionEnd + 1;
        end
    end
    TrialDataUtilities.Plotting.hideInLegend(hs);
end

function [ha] = shadeSimple(axh, x, y1, y2, z, varargin)

p = inputParser();
p.addParameter('FaceColor', [0.8 0.8 1], @(x) true);
p.addParameter('EdgeColor', 'none', @(x) true);
p.addParameter('alpha', 1, @isscalar);
p.addParameter('Clipping', true, @islogical);
p.KeepUnmatched = false;
p.parse(varargin{:});

faceColor = Neuropixel.Utils.convertColorsToMatrix(p.Results.FaceColor);
edgeColor = Neuropixel.Utils.convertColorsToMatrix(p.Results.EdgeColor);

xv = [x, fliplr(x)];
yv = [y1, fliplr(y2)];
zv = z * ones(size(xv));

ha = patch(xv, yv, zv, 'k', 'Parent', axh);
set(ha, 'FaceColor', faceColor, ...
    'EdgeColor', edgeColor, 'Parent', axh, 'FaceAlpha', p.Results.alpha, 'Clipping', p.Results.Clipping);

% hide shading from legend
set(get(get(ha, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

% set(axh, 'Layer', 'top')

end
