function data = demeanOverTime(imec, data, chIdx, timeIdx) %#ok<INUSD>
    ch_mask = ismember(chIdx, imec.connectedChannels);
    data(ch_mask, :) = data(ch_mask, :) - median(data(ch_mask),2); % make median of each channel 0 over tim
end
