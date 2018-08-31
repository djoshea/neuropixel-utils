function showFirstInLegend(h, name)
% showInLegend(h, names)
% show first object in h in default legend with specified name
% use legend(axh, 'show') to activate default legend

    if isempty(name)
        hideInLegend(h);
    else
        h = h(isGraphicsHandle(h));
        if ~isempty(h)
            showInLegend(h(1), name);
            if numel(h) > 1
                hideInLegend(h(2:end));
            end
        end
    end
end

function showInLegend(h, names)
% showInLegend(h, names)
% show object h in legend with specified by default
% use legend(axh, 'show') to activate default legend
% names is either char (for scalar h), or cellstr

    if nargin >= 2 && ~isempty(names) && ischar(names)
        names = repmat({names}, numel(h), 1);
    end

    for i = 1:numel(h)
        if ~isGraphicsHandle(h(i)), continue; end
        ann = get(h(i), 'Annotation');
        leg = get(ann, 'LegendInformation');
        set(leg, 'IconDisplayStyle', 'on');
        
        if ~isempty(names)
            set(h(i), 'DisplayName', names{i});
        end
    end

end

function hideInLegend(h)
    % prevent object h from appearing in legend by default
    for i = 1:numel(h)
        if ~isGraphicsHandle(h(i)), continue; end
        ann = get(h(i), 'Annotation');
        leg = get(ann, 'LegendInformation');
        set(leg, 'IconDisplayStyle', 'off');
    end
end

function mask = isGraphicsHandle(h)
    if isobject(h)
        mask = isgraphics(h);
    else
        mask = ~isnan(h);
    end
end