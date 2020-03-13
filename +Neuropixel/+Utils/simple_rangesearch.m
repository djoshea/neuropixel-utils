function idx = simple_rangesearch(sorted_list, lookups, window)
% simple alternative to rangesearch where sorted_list is a 1d sorted list. Finds indices of sorted_list that are within lookups(i) + window(1) : lookups(i) + window(2)

    N = numel(sorted_list);
    sorted_list = cast(sorted_list, 'like', window);
    
    [sorted_lookups, idx_sort_lookups] = sort(lookups);
    K = numel(sorted_lookups);

    spotlight_left = 1;
    c_idx = cell(K, 1);
    for k = 1:K
        lookup = cast(sorted_lookups(k), class(window));
        win_start = lookup + window(1);
        win_stop = lookup + window(2);

        [c_idx{k}, spotlight_left] = find_window_bounds(win_start, win_stop, spotlight_left);
    end
    
    % undo sorting of lookups
    idx(idx_sort_lookups) = c_idx;

    function [idx, start_at] = find_window_bounds(win_start, win_stop, start_at)
        % advance to start of window
        for iV = start_at:N
            if sorted_list(iV) >= win_start
                break;
            end
        end
        start_at = iV;

        idx_last = start_at - 1;
        % find end of window
        for iV = start_at:N
            if sorted_list(iV) <= win_stop
                idx_last = iV;
            else
                break;
            end
        end

        idx = start_at:idx_last;

    end

end