function [out, transformInfo] = plotStackedTraces(tvec, data_txcxl, varargin)
% data is time x channels x layers

p = inputParser();
p.addParameter('style', 'traces', @(x) ischar(x) || isstring(x)); 
p.addParameter('transformInfo', [], @(x) isempty(x) || isstruct(x));

p.addParameter('colors', [], @(x) true); % over channels
p.addParameter('layerColors', [], @(x) true); % over layers (specify one or the other)
p.addParameter('normalizeEach', false, @islogical);
p.addParameter('normalizeMask', [], @isvector);
p.addParameter('quantile', 1, @isscalar);
p.addParameter('lineArgs', {}, @iscell);
p.addParameter('LineWidth', 1, @isscalar);
p.addParameter('LineOpacity', 1, @isscalar);
p.addParameter('RidgeOpacity', 1, @isscalar);
p.addParameter('gain', 1, @isscalar);
p.addParameter('labels', [], @(x) true);
p.addParameter('channel_ids', [], @(x) true); % used for data tips, c x l or c x 1
p.addParameter('invertChannels', false, @islogical);

% common data tip for all traces
% value is of size(data) (t x c x l) and singleton dimensions will be automatically expanded
p.addParameter('dataTipLabel', 'overlay', @(x) isempty(x) || ischar(x) || isstring(x));
p.addParameter('dataTipValues', [], @(x) true); 
p.addParameter('dataTipFormat', '', @ischar);

p.addParameter('showChannelDataTips', true, @islogical);
p.addParameter('showOverlayDataTips', true, @islogical);

p.addParameter('colorOverlayLabels', [], @(x) isempty(x) || isnumeric(x)); % a label matrix indexing into colorOverlapMap of size data_txc
p.addParameter('colorOverlayMap', [], @(x) isempty(x) || ismatrix(x)); % a Nx3 colormap
p.addParameter('colorOverlayDataTipLabel', 'overlay', @(x) isempty(x) || ischar(x) || isstring(x));
p.addParameter('colorOverlayDataTipValues', [], @(x) true); % nlabels x 1
p.addParameter('colorOverlayDataTipFormat', '', @ischar);
p.parse(varargin{:});

showChannelDataTips = p.Results.showChannelDataTips;
showOverlayDataTips = p.Results.showOverlayDataTips;

style = string(p.Results.style);
assert(ismember(style, ["traces", "ridgeline", "heatmap"]));

data = data_txcxl;
if isinteger(data)
    data = single(data);
end
nTraces = size(data, 2);
nLayers = size(data, 3);

traceColors = p.Results.colors;
if ~isempty(traceColors) && size(traceColors, 1) == 1
    traceColors = repmat(traceColors, nTraces, 1);
end
layerColors = p.Results.layerColors; 
if ~isempty(layerColors) && size(layerColors, 1) == 1
    layerColors = repmat(layerColors, nLayers, 1);
end

if isempty(traceColors) && isempty(layerColors)
    if nLayers > 1
        %layerColors = Neuropixel.Utils.phy_cluster_colors(nLayers);
        layerColors = Neuropixel.Utils.linspecer(nLayers);
    else
        traceColors = repmat([0 0 0], nTraces, 1);
    end
end

tform = p.Results.transformInfo;
if isempty(tform)
    tform = struct();
end

% optionally keep a subset of channels intact
if ~isfield(tform, 'cmask')
    if ~isempty(p.Results.normalizeMask)
        tform.cmask = p.Results.normalizeMask;
    else
        tform.cmask = true(size(data, 2), 1);
    end
end
cmask = tform.cmask;

if ~isfield(tform, 'bias')
    tform.bias = - mean(data(:, cmask, :), 1, 'omitnan');
end
data(:, cmask, :) = data(:, cmask, :) + tform.bias;

if ~isfield(tform, 'normBy')
    if p.Results.normalizeEach
        tform.normBy = quantile(abs(data(:, cmask, :)), p.Results.quantile, 1);
    else
        vals = data(:, cmask, :);
        tform.normBy = quantile(abs(vals(:)), p.Results.quantile);
    end
end

if isscalar(tform.normBy) && tform.normBy ~= 0
    data(:, cmask, :) = data(:, cmask, :) ./ tform.normBy;
elseif numel(tform.normBy) > 1
    data(:, cmask & tform.normBy ~= 0, :) = data(:, cmask & tform.normBy ~= 0, :) ./ tform.normBy;
end

if ~isfield(tform, 'gain')
    tform.gain = p.Results.gain; 
end
data(:, cmask, :) = data(:, cmask, :) * tform.gain;
data(:, ~cmask, :) = data(:, ~cmask, :) * 0.95; % this keeps logical bits from touching

multipliers = nan(size(data, 2), 1);
multipliers(cmask) = tform.normBy * tform.gain;

if style == "traces" || style == "ridgeline"
    if ~isfield(tform, 'offsets')
        if p.Results.invertChannels % first is at bottom
            tform.offsets = 0:nTraces-1;
        else
            tform.offsets = nTraces-1:-1:0; % default, first is at top
        end
    end
else
    tform.offsets = zeros(1, nTraces);
end

data = data + tform.offsets; % channels along sercond dimension

washolding = ishold;
labels = string(p.Results.labels);

% possibly different channel_ids for each layer
channel_ids = p.Results.channel_ids;
if ~isempty(channel_ids)
    if size(channel_ids, 2) == 1
        channel_ids = repmat(channel_ids, 1, nLayers);
    end
    assert(size(channel_ids, 1) == nTraces);
end

hvec = gobjects(nTraces, nLayers);
dataTipValue = p.Results.dataTipValues;
if ~isempty(dataTipValue)
    if size(dataTipValue, 1) == 1
        dataTipValue = repmat(dataTipValue, size(data, 1), 1, 1);
    end
    if size(dataTipValue, 2) == 1
        dataTipValue = repmat(dataTipValue, 1, size(data, 2), 1);
    end
    if size(dataTipValue, 3) == 1
        dataTipValue = repmat(dataTipValue, 1, 1, size(data, 3));
    end
    if isempty(p.Results.dataTipFormat)
        dtipfmatArg = {};
    else
        dtipfmatArg = {p.Results.dataTipFormat};
    end
end

switch style
    case {"traces", "ridgeline"}
        for iR = 1:nTraces
            dmat = squeeze(data(:, iR, :));
            if ~any(~isnan(dmat), 'all'), continue, end
            
            if style == "traces"
                hvec(iR, :) = plot(tvec, dmat, '-', 'Color', [0 0 0 p.Results.LineOpacity], 'LineWidth', p.Results.LineWidth, p.Results.lineArgs{:});
                if ~isempty(traceColors)
                    % color by trace
                    set(hvec(iR, :), 'Color', [traceColors(iR, :) p.Results.LineOpacity]);
                else
                    % color by layers
                    for iL = 1:nLayers
                        hvec(iR, iL).Color = [layerColors(iL, :) p.Results.LineOpacity];
                    end
                end
            else
                hvec(iR, :) = area(tvec, dmat, tform.offsets(iR), 'ShowBaseLine', 'off', 'FaceColor', 'w', 'EdgeAlpha', p.Results.LineOpacity, 'FaceAlpha', p.Results.RidgeOpacity);
                if ~isempty(traceColors)
                    % color by trace
                    set(hvec(iR, :), 'EdgeColor', traceColors(iR, :));
                else
                    % color by layers
                    for iL = 1:nLayers
                        hvec(iR, iL).EdgeColor = layerColors(iL, :);
                    end
                end
            end
            

            if ~verLessThan('matlab', '9.6.0') 
                if ~isempty(channel_ids) && showChannelDataTips % R2019a 
                    for iL = 1:nLayers
                        nPoints = numel(hvec(iR, iL).XData);
                        hvec(iR, iL).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('channel', repmat(double(channel_ids(iR, iL)), nPoints, 1), '%d');
                    end
                end

                if ~isempty(dataTipValue)
                    for iL = 1:nLayers
                        hvec(iR, iL).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow(p.Results.dataTipLabel, double(dataTipValue(:, iR, iL)), dtipfmatArg{:});
                    end
                end
            end

            hold on;

        %     if ~isempty(labels)
        %         text(tvec(1) - (tvec(2) - tvec(1)), offsets(iR), labels{iR}, 'Background', 'none', ...
        %             'HorizontalAlign', 'right', 'VerticalAlign', 'bottom', 'Color', colors(iR, :));
        %     end
        end

    case 'heatmap'
        assert(nLayers == 1);
        Neuropixel.Utils.pmatbal(data');
end

% do color overlays
if ~isempty(p.Results.colorOverlayLabels)
    labelMat = p.Results.colorOverlayLabels;
    
    dataTipLabel = p.Results.colorOverlayDataTipLabel;
    dataTipValuesByLabel = p.Results.colorOverlayDataTipValues;
    dataTipFormat = p.Results.colorOverlayDataTipFormat;
    
    assert(isequal(size(labelMat), size(data)));
    maxLabel = max(labelMat(:));
    
    labelCmap = p.Results.colorOverlayMap;
    if isempty(labelCmap)
        labelCmap = Neuropixel.Utils.distinguishable_colors(maxLabel, {'k', 'w'});
    end
    
    hLabel = cell(maxLabel, nLayers);
    for iU = 1:maxLabel
        mask = labelMat == iU;
        if ~any(mask(:)), continue, end
        dataMasked = data;
        dataMasked(~mask) = NaN;
        
        for iL = 1:nLayers
            t_mask = any(labelMat(:, :, iL) == iU, 2);
            % dilate this time mask by one to the right, to preserve the NaNs in between valid regions
            t_mask = imdilate(t_mask, [1; 1]);
            ch_mask = any(labelMat(t_mask, :, iL) == iU, 1);
            
            hLabel{iU, iL} = plot(tvec(t_mask), dataMasked(t_mask, ch_mask, iL), '-', 'Color', labelCmap(iU, :), 'LineWidth', p.Results.LineWidth, p.Results.lineArgs{:});
            
            if ~verLessThan('matlab', '9.6.0') && showOverlayDataTips && ~isempty(dataTipValuesByLabel) % R2019a 
                for iH = 1:numel(hLabel{iU, iL})
                    nPoints = numel(hLabel{iU, iL}(iH).XData);
                    label_this = repmat(dataTipValuesByLabel(iU), nPoints, 1);
                    if isnumeric(label_this)
                        label_this = double(label_this);
                    end
                    if isempty(dataTipFormat)
                        fmatArg = {};
                    else
                        fmatArg = {dataTipFormat};
                    end
                    
                    hLabel{iU, iL}(iH).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow(dataTipLabel, label_this, fmatArg{:});
                    channel_ids_this = channel_ids(ch_mask, iL);
                    hLabel{iU, iL}(iH).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('channel', repmat(double(channel_ids_this(iH)), nPoints, 1), '%d');
                end
            end
        end
    end
    
else
    hLabel = {};
end

if ~isempty(labels)
    [sorted_offsets, sort_idx] = sort(tform.offsets);
    set(gca, 'YTick', sorted_offsets, 'YTickLabels', labels(sort_idx));
end


out.tvec = tvec;
out.data = data;
out.hLines = hvec;
out.multipliers = multipliers;
out.hLabel = hLabel;

transformInfo = tform;

axis tight;
box off;

if ~isempty(traceColors)
    cmap = traceColors;
else
    cmap =layerColors;
end

ax = gca;
ax.TickDir = 'out';
if style == "traces" || style == "ridgeline"
    ax.ColorSpace.Colormap = cmap;
    if size(cmap, 1) == 1
        ax.CLim = [0.9 1.1];
    else
        ax.CLim = [1 size(cmap, 1)];
    end
end

if ~washolding
    hold off;
end


end