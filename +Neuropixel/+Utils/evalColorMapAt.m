function ceval = evalColorMapAt(cmap, at, clim)
    if isa(cmap, 'function_handle')
        cmap = cmap(1000);
    end
    if nargin < 3
        clim = [0 1];
    end
    
    clim = cast(clim, 'like', cmap);
    if clim(2) == clim(1)
        clim(2) = clim(2) + 0.1;
    end
    
    N = size(cmap, 1);
    at = cast(at, 'like', cmap);
    
    cmapOrigAt = (0:N-1) / (N-1) * (clim(2)-clim(1)) + clim(1);
    ceval = interp1(cmapOrigAt, cmap, at);
    ceval = max(min(ceval, 1), 0);
end