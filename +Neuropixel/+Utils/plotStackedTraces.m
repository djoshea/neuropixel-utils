function plotStackedTraces(tvec, data_txc, varargin)

p = inputParser();
p.addParameter('colors', [], @(x) true);
p.addParameter('normalizeEach', false, @isvarargin);
p.addParameter('lineArgs', {}, @iscell);
p.addParameter('height', 0.95, @isscalar);
p.addParameter('labels', [], @(x) true);
p.parse(varargin{:});

data = data_txc;
nTraces = size(data, 2);

colors = p.Results.colors;
if isempty(colors)
    colors = zeros(nTraces, 3);
end

data = data - min(data, [], 1, 'omitnan');

if p.Results.normalizeEach
    data = data ./ max(data,[], 1, 'omitnan');
else
    data = data ./ max(data(:), [], 'omitnan');
end
data = data * p.Results.height;

offsets = nTraces-1:-1:0;
data = data + offsets;

washolding = ishold;
labels = string(p.Results.labels);

for iR = 1:nTraces
    plot(tvec, data(:, iR), '-', 'Color', colors(iR, :), p.Results.lineArgs{:});
    hold on;
    
    if ~isempty(labels)
        text(tvec(1), offsets(iR), labels{iR}, 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'Color', colors(iR, :));
    end
end

axis tight;
box off;

if ~washolding
    hold off;
end

end