function data = commonAverageReferenceEachBank(imec, data, chIds, timeIdx, extraArg) %#ok<INUSL>
    if nargin < 5
        extraArg = struct();
    end
    
    chanMaskGood = ismember(chIds, imec.goodChannels);

    % subtract median of each channel over time
    data(chanMaskGood, :) = bsxfun(@minus, data(chanMaskGood, :), median(data(chanMaskGood, :), 2)); 

    % subtract median across good channels
    data(chanMaskGood, :) = bsxfun(@minus, data(chanMaskGood, :), median(data(chanMaskGood, :), 1));
    
    % then do Siegle et al. 2019 style median subtraction over simultaneously acquired samples
    if contains(imec.channelMap.name, 'phase3a', 'IgnoreCase', true)
        for n = 1:24
            chIdsThisGroup = n:24:384;
            chanMaskThisGroup = ismember(chIds, chIdsThisGroup) & chanMaskGood;
            data(chanMaskThisGroup, :) = data(chanMaskThisGroup, :) - median(data(chanMaskThisGroup, :), 1, 'omitnan');
        end
    end
    
    % optionally do high pass filter on data as vanilla KS 2 would, this is useful to eliminate post-stim baseline drift
    if isfield(extraArg, 'hp_filter') && extraArg.hp_filter
        [b, a] = butter(extraArg.hp_filter_half_order, extraArg.hp_filter_corner/extraArg.fs*2, 'high');
        data(chanMaskGood, :) = filter(b, a, data(chanMaskGood, :), [], 2); % causal forward filter
        data(chanMaskGood, :) = fliplr(filter(b, a, fliplr(data(chanMaskGood, :)), [], 2)); % acausal reverse filter
    end
end

