function cmap = phy_cluster_colors(N)
% I don't know what this is, but it looks graphics-related
%
% cmap = npxutils.internal.graphics.phy_cluster_colors(N)

% copied from phy.utils._color.ColorSelector
cmap = [
    8, 146, 252
    255, 2, 2
    240, 253, 2
    228, 31, 228
    2, 217, 2
    255, 147, 2
    
    212, 150, 70
    205, 131, 201
    201, 172, 36
    150, 179, 62
    95, 188, 122
    129, 173, 190
    231, 107, 119
    ] ...
    ./ 255;

if nargin > 0
    nmap = size(cmap, 1);
    if nmap < N
        cmap = repmat(cmap, ceil(N / nmap), 1);
    end
    cmap = cmap(1:N, :);
end

end