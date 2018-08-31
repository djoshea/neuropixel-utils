function data = demeanOverTime(imec, data)
    ch_mask = imec.connectedChannels;
    data(ch_mask, :) = data(ch_mask, :) - median(data(ch_mask),2); % make median of each channel 0 over tim
end
