function image_pos = computeWaveformImageCentroid(image_data, channel_inds_by_image, ch_pos, varargin)
% image_data is nImages x nTime x nChannelsPerImage
% channel_inds_by_image is nImages x nChannelsPerImage lookup of inds into ch_pos (in [1 nChannels])
% ch_pos is nChannels x n_spatial_dims (typically 2 or 3)
%
% image_pos is nImages x n_spatial_dims
%
% methods:
%   - com_via_range : takes the weighted centroid using the range (max-min) on each channel, the original technique I implemented
p = inputParser();
p.addParameter('method', 'com_via_range', @ischar);
p.addParameter('relativeThreshold', 0.1, @isscalar);
p.parse(varargin{:});

method = p.Results.method;
relThresh = p.Results.relativeThreshold;

[nImages, nTime, nChPerImage] = size(image_data);
if size(channel_inds_by_image, 1) == 1
    channel_inds_by_image = repmat(channel_inds_by_image, nImages, 1);
end
assert(size(channel_inds_by_image, 1) == nImages);
assert(size(channel_inds_by_image, 2) == nChPerImage);

nSpatialDims = size(ch_pos, 2);

% reshape (nImages x nChannelsPerImage) x (spatial_dim) --> nImages x nChannelsPerImage x spatial_dim
image_ch_pos = reshape(ch_pos(channel_inds_by_image(:), :), [nImages, nChPerImage, nSpatialDims]);

switch method
    case 'com_via_range'
        % compute the center of mass taking the amplitude as the range of each channel individually
        amp_data_raw = Neuropixel.Utils.TensorUtils.squeezeDims(max(image_data,[],2) - min(image_data,[],2), 2); % nImages x nChannelsPerImage
        peak_amp = max(amp_data_raw, [], 2);
        amp_data_thresh = amp_data_raw;
        amp_data_thresh(amp_data_raw < peak_amp * relThresh) = 0;
        
        % weighted sum over channels in each dimension
        % amp (nImages x nChPerImage) .* pos (nImages x nChPerImage x spatialDim) over ch, divide by sum amps over ch
        % --> nImages x 1 x spatialDim --> nIamges x spatialDim
        image_pos = Neuropixel.Utils.TensorUtils.squeezeDims(sum(amp_data_thresh .* image_ch_pos, 2) ./ ...
            sum(amp_data_thresh, 2), 2);

    case 'com_via_provided_amplitudes'
        % compute the center of mass taking the amplitude as given by the single timepoint in image_data
        assert(nTime == 1);
        amp_data_raw = Neuropixel.Utils.TensorUtils.squeezeDims(image_data, 2);
        peak_amp = max(amp_data_raw, [], 2);
        
        amp_data_thresh = amp_data_raw;
        amp_data_thresh(image_data < peak_amp * relThresh) = 0;

        image_pos = Neuropixel.Utils.TensorUtils.squeezeDims(sum(amp_data_thresh .* image_ch_pos, 2) ./ ...
            sum(amp_data_thresh, 2), 2);
        
    case 'com_best_timepoint'
        % % compute the center of mass taking the amplitude as the range of each channel individually
        [~, best_linear_ind] = max(abs(image_data(:, :)), [], 2);
        [best_time, best_ch] = ind2sub([nTime, nChPerImage], best_linear_ind);
        
        amp_data_raw = nan(nImages, nChPerImage, 'like', image_data);
        for iI = 1:nImages
            % take values at that timepoint with the same sign as the extremum
            val_best = image_data(iI, best_time(iI), best_ch(iI));
            if val_best < 0
                amp_data_raw(iI, :) = max(-image_data(iI, best_time(iI), :), 0);
            else
                amp_data_raw(iI, :) = max(image_data(iI, best_time(iI), :), 0);
            end
        end
        
        peak_amp = max(amp_data_raw, [], 2);
        amp_data_thresh = amp_data_raw;
        amp_data_thresh(amp_data_raw < peak_amp * relThresh) = 0;
        
        % weighted sum over channels in each dimension
        % amp (nImages x nChPerImage) .* pos (nImages x nChPerImage x spatialDim) over ch, divide by sum amps over ch
        % --> nImages x 1 x spatialDim --> nIamges x spatialDim
        image_pos = Neuropixel.Utils.TensorUtils.squeezeDims(sum(amp_data_thresh .* image_ch_pos, 2) ./ ...
            sum(amp_data_thresh, 2), 2);

    otherwise
        error('Unknown centroid method %s', method);
        
end


end