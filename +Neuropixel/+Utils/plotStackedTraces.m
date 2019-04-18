function plotStackedTraces(tvec, data_txc, varargin)

p = inputParser();
p.addParameter('colors', [], @(x) true);
p.addParameter('normalizeEach', false, @islogical);
p.addParameter('normalizeMask', [], @isvector);
p.addParameter('quantile', 1, @isscalar);
p.addParameter('lineArgs', {}, @iscell);
p.addParameter('gain', 10, @isscalar);
p.addParameter('labels', [], @(x) true);
p.parse(varargin{:});

data = data_txc;
nTraces = size(data, 2);

colors = p.Results.colors;
if isempty(colors)
    colors = zeros(nTraces, 3);
end

% optionally keep a subset of channels intact
cmask = p.Results.normalizeMask;
if isempty(cmask)
    cmask = true(size(data, 2), 1);
end
data(:, cmask) = data(:, cmask) - mean(data(:, cmask), 1, 'omitnan');


if p.Results.normalizeEach
    normBy = quantile(abs(data(:, cmask)), p.Results.quantile, 1);
    data(:, cmask) = data(:, cmask) ./ normBy;
else
    vals = data(:, cmask);
    normBy = quantile(abs(vals(:)), p.Results.quantile);
    data(:, cmask) = data(:, cmask) ./ normBy;
end
data(:, cmask) = data(:, cmask) * p.Results.gain;

offsets = nTraces-1:-1:0;
data = data + offsets;

washolding = ishold;
labels = string(p.Results.labels);

for iR = 1:nTraces
    plot(tvec, data(:, iR), '-', 'Color', colors(iR, :), p.Results.lineArgs{:});
    hold on;
    
%     if ~isempty(labels)
%         text(tvec(1) - (tvec(2) - tvec(1)), offsets(iR), labels{iR}, 'Background', 'none', ...
%             'HorizontalAlign', 'right', 'VerticalAlign', 'bottom', 'Color', colors(iR, :));
%     end
end

if ~isempty(labels)
    set(gca, 'YTick', fliplr(offsets), 'YTickLabels', flipud(labels));
end

axis tight;
box off;
set(gca, 'TickDir', 'out');

if ~washolding
    hold off;
end


end