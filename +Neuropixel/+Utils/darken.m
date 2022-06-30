function cmap = darken(cmap, amount)
    if nargin < 2
        amount = 1;
    end
    
    Kn = 18;
    lab = Neuropixel.Utils.convert('RGB->Lab', cmap);

    lab(:, 1) = lab(:, 1) - Kn * amount;

    cmap = Neuropixel.Utils.convert('Lab->RGB', lab);

end