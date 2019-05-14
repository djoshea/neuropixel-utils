function [h, hcbar] = pmatbal(m, varargin)
% visualize a matrix using pcolor but with a symmetric red white blue colormap

[h, hcbar] = Neuropixel.Utils.pmat(m, varargin{:});
Neuropixel.Utils.cmocean('balance');
L = max(abs(m(:)));
caxis([-L L]);
