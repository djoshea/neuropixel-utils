function data = commonAverageReferenceEachBank(imec, data, chIds, timeIdx) %#ok<INUSD>
    chanMask = ismember(chIds, imec.goodChannels);

    % subtract median of each channel over time
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 2)); 

    % subtract median across good channels
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 1));
    
    % then do Siegle et al. 2019 style median subtraction over simultaneously acquired samples
    if contains(imec.channelMap.name, 'phase3a', 'IgnoreCase', true)
        for n = 1:24
            chIdsThisGroup = n:24:384;
            chanMaskThisGroup = ismember(chIds, chIdsThisGroup) & chanMask;
            data(chanMaskThisGroup, :) = data(chanMaskThisGroup, :) - median(data(chanMaskThisGroup, :), 1, 'omitnan');
        end
    end
end

