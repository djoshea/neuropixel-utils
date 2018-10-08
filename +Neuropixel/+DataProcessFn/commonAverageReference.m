function data = commonAverageReference(imec, data, chIdx, timeIdx) %#ok<INUSD>
    chanMask = ismember(chIdx, imec.goodChannels);

    % subtract median of each channel over time
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 2)); 

    % subtract median across good channels
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 1));
end

