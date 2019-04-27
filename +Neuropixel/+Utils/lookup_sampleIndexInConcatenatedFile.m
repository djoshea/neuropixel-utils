function [fileInds, origSampleInds] = lookup_sampleIndexInConcatenatedFile(starts, inds)
   % convert from a time index in this file to which of the concatenated files that index came from

   [fileInds, origSampleInds] = deal(nan(size(inds)));
   stops = [Neuropixel.Utils.makecol(starts(2:end)); Inf];
   
   for i = 1:numel(starts)
      mask = inds >= starts(i) & inds < stops(i);
      fileInds(mask) = i;
      origSampleInds(mask) = inds(mask) - starts(i);
   end
   
end