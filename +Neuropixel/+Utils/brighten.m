function cmap = brighten(cmap, amount)
    if nargin < 2
        amount = 1;
    end
    
    cmap = Neuropixel.Utils.darken(cmap, -amount);
end