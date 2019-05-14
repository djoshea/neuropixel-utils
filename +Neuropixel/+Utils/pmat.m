function [h, hcbar] = pmat(m, varargin)
% visualize a matrix using pcolor

p = inputParser();
p.addParameter('x', [], @(x) isvector(x));
p.addParameter('y', [], @(x) isvector(x));
p.addParameter('addColorbar', true, @islogical);
p.parse(varargin{:});

cla;

m = squeeze(m);

if isvector(m)
    m = repmat(makerow(m), 2, 1);
end

if islogical(m)
    m = double(m);
end

if ndims(m) > 2 %#ok<ISMAT>
    warning('Selecting (:, :, 1) of tensor to display');
    m = m(:, :, 1);
end

% add an extra row onto m
addRowCol = @(v) [v, v(:, end)+diff(v(:, end-1:end), 1, 2); ...
    v(end, :) + diff(v(end-1:end, :), 1, 1), 2*v(end, end)-v(end-1, end-1)];

if isempty(p.Results.x)
    x = 0.5:size(m, 2)-0.5;
else
    x = p.Results.x;
    dx = diff(x);
    x = x - dx([1:end end])/2;
end
if isempty(p.Results.y)
    y = 0.5:size(m, 1)-0.5;
else
    y = p.Results.y;
    dy = diff(y);
    y = y - dy([1:end end])/2;
end

[X, Y] = meshgrid(x, y);
        
% need an extra row and column because of the way that pcolor works
m = addRowCol(m);
X = addRowCol(X);
Y = addRowCol(Y);

h = pcolor(X,Y, m);

set(h, 'EdgeColor', 'none');
Neuropixel.Utils.cmocean('haline');

if p.Results.addColorbar
    hcbar = colorbar;
    box(hcbar, 'off');
    set(hcbar, 'TickLength', 0);
else
    hcbar = [];
end

box off
axis ij
axis on;

axis on;
set(gca, 'TickLength', [0 0], 'XAxisLocation', 'top');
axis tight;
box on;

