function which_window = discretize_windows(values, windows)
% values is N x 1, windows is K x 2 with [lower, upper] bound on each window
% which_window is N x 1 indices in range 1:K or NaN if not in window
% windows must be non-overlapping, but need not be in sorted order
%
% by default, the upper and lower edges of a window are included in a window.
% but if two windows share a common edge, the value will be assigned to the later window (whose lower bound matches the value)
%
% test:
% >> Neuropixel.Utils.discretize_windows([0.5 1 1.5 2 2.5 3 3.5 4.5], [0.5 1.5; 1.5 2.5; 3.5 4.5])
%
% ans =
%
%     1     1     2     2     2   NaN     3     3
%

    % sort windows
    nWindows = size(windows, 1);
    assert(size(windows, 2) == 2);
    [sorted_windows, idx_sort_windows] = sortrows(windows, 1);

    % extract monotonically increasing bin edges with windows and gaps between windows side by side
    edges = sorted_windows';
    edges = double(edges(:));

    % add small epsilon to the right edge of each window so that each bin includes its right edge, provided
    % that it has a gap between it and the next window
    window_gaps = diff(edges);
    window_gaps = window_gaps(2:2:end);
    inds_add_eps = 2:2:2*nWindows;
    inds_add_eps(window_gaps < eps) = [];
    edges(inds_add_eps) = edges(inds_add_eps) + eps(edges(inds_add_eps));

    % discretize into bins
    temp_bins = discretize(values, edges, 'IncludedEdge', 'left');

    % remove even values which lie between bins
    mask_in_window = mod(temp_bins, 2) == 1;
    sorted_bin = (temp_bins+1) / 2;
    which_window = nan(size(values));
    which_window(mask_in_window) = idx_sort_windows(sorted_bin(mask_in_window));

end