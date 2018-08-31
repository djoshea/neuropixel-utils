function data = commonAverageReference(imec, data)
    chanMask = ismember(1:imec.nChannelsMapped, imec.goodChannels);

    data = bsxfun(@minus, data, median(data,2)); % make median of each channel 0

    % subtract median of good channels
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 1));
end

