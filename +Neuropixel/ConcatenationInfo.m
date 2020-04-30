classdef ConcatenationInfo < handle
    properties
        fs
        samplesPreShift (:, 1) uint64
        timeShifts (:, 1) % Neuropixel.TimeShiftSpec
        names (:, 1) string
        gains (:, 1) 
        multipliers (:, 1) 
        adcBits (:, 1)
        ranges (:, 2)
        nSamples
    end
    
    properties(Dependent)
        nDatasets
        samples % nDatasets x 1 number of samples in resultant file, factoring in time shifts
        startIdx
        stopIdx
        
        scaleToUvs
    end
    
    methods
        function ci = ConcatenationInfo(imec, mode, meta)
            if nargin > 0
                if ~exist('meta', 'var'), meta = []; end
                
                switch mode
                    case 'ap'
                        if isempty(meta), meta = imec.readAPMeta(); end
                        ci.fs = imec.fsAP;
                        ci.nSamples = imec.nSamplesAP;
                        ci.samplesPreShift = getor(meta, 'concatenatedSamples', imec.nSamplesAP);
                    case 'lf'
                        if isempty(meta), meta = imec.readLFMeta(); end
                        ci.fs = imec.fsLF;
                        ci.nSamples = imec.nSamplesLF;
                        ci.samplesPreShift = getor(meta, 'concatenatedSamples', imec.nSamplesLF);
                    otherwise
                        error ('unknown mode');
                end
                
                nC = numel(ci.samplesPreShift);
                
                timeShifts = getor(meta, 'concatenatedTimeShifts', []);
                if ~isempty(timeShifts)
                    if isnumeric(timeShifts)
                        % single time shift, deserialized as single 1 x 3 numeric matrix
                        ci.timeShifts = Neuropixel.TimeShiftSpec.from_matrix(timeShifts, ci.samplesPreShift);
                    else
                        % expecting mat2str serialized [1 2 3;4 5 6] type arrays, joined by ;
                        matches = string(regexp(timeShifts, '\s*(\[[\d ;,]*])\s*;?\s*', 'tokens'))';
                        ci.timeShifts = arrayfun(@(str, nSamples) Neuropixel.TimeShiftSpec.from_string(str, nSamples), matches, ci.samplesPreShift);
                    end
                    assert(sum([ci.timeShifts.nShiftedTimes]) == ci.nSamples, 'Mismatch between time shifts and total number of samples');
                end
                
                namescat = getor(meta, 'concatenated', []);
                if ~isempty(namescat)
                    if contains(namescat, ';')
                        ci.names = string(Neuropixel.Utils.makecol(strsplit(namescat, ';')));
                    else
                        ci.names = string(Neuropixel.Utils.makecol(strsplit(namescat, ':')));
                    end
                else
                    ci.names = strings(nC, 1);
                end
                
                ci.gains = getor(meta, 'concatenatedGains', imec.apGain * ones(nC, 1));
                ci.multipliers = getor(meta, 'concatenatedMultipliers', ones(nC, 1));
                ci.adcBits = getor(meta, 'concatenatedAdcBits', imec.adcBits * ones(nC, 1));
                
                switch mode
                    case 'ap'
                        rmin = Neuropixel.Utils.makecol(getor(meta, 'concatenatedAiRangeMin', imec.apRange(1) * ones(nC, 1)));
                        rmax = Neuropixel.Utils.makecol(getor(meta, 'concatenatedAiRangeMax', imec.apRange(2) * ones(nC, 1)));
                    case 'lf'
                        rmin = Neuropixel.Utils.makecol(getor(meta, 'concatenatedAiRangeMin', imec.lfRange(1) * ones(nC, 1)));
                        rmax = Neuropixel.Utils.makecol(getor(meta, 'concatenatedAiRangeMax', imec.lfRange(2) * ones(nC, 1)));
                    otherwise
                        error ('unknown mode');
                end 
                ci.ranges = cat(2, rmin, rmax);
            end
            
            function v = getor(s, fld, def)
                if isfield(s, fld)
                    v = s.(fld);
                else
                    v = def;
                end
            end
        end
    end
    
    methods(Static)
        function ci = inferFromOtherBand(imec, bandInfer, bandRef)
            ci = Neuropixel.ConcatenationInfo();
            
            bandInfer = string(bandInfer);
            bandRef = string(bandRef);
            switch bandInfer
                case "ap"
                    assert(imec.hasSourceAP, 'Must have source AP datasets to infer AP concatenationInfo');
                    fsInfer = imec.sourceDatasets(1).fsAP;
                    ci.samplesPreShift = arrayfun(@(src) src.nSamplesAP, imec.sourceDatasets);
                    
                case "lf"
                    assert(imec.hasSourceLF, 'Must have source AP datasets to infer LF concatenationInfo');
                    fsInfer = imec.sourceDatasets(1).fsLF;
                    ci.samplesPreShift = arrayfun(@(src) src.nSamplesLF, imec.sourceDatasets);
                    
                otherwise
                    error('Unknown bandInfer %s', bandInfer);
            end
            
            assert(bandInfer ~= bandRef, 'bandInfer and bandRef must be different');
            
            switch bandRef
                case "ap"
                    assert(imec.hasAP, 'ImecDataset lacks AP band');
                    ciRef = imec.concatenationInfoAP;
                    fsRef = imec.fsAP;
                case "lf"
                    assert(imec.hasLF, 'ImecDataset lacks LF band');
                    ciRef = imec.concatenationInfoLF;
                    fsRef = imec.fsLF;
                otherwise
                    error('Unknown bandSource %s', bandSource);
            end
            
            function samplesInfer = samplesRefToInfer(samplesRef)
                samplesInfer = cast(floor(double(samplesRef) ./ double(fsRef) .* double(fsInfer)), 'like', samplesRef);
            end
            
            ci.fs = fsInfer;
            if ~isempty(ciRef.timeShifts)
                ci.timeShifts = copy(ciRef.timeShifts);
            else
                ci.timeShifts = ciRef.timeShifts;
            end
            nSamplesInfer = uint64(0);
            for iTS = 1:numel(ci.timeShifts)
                ci.timeShifts(iTS).idxStart = samplesRefToInfer(ci.timeShifts(iTS).idxStart);
                ci.timeShifts(iTS).idxStop = samplesRefToInfer(ci.timeShifts(iTS).idxStop);
                ci.timeShifts(iTS).idxShiftStart = samplesRefToInfer(ci.timeShifts(iTS).idxShiftStart);
                
                nSamplesInfer = nSamplesInfer + ci.timeShifts(iTS).nShiftedTimes;
            end
            ci.nSamples = nSamplesInfer;
            ci.names = ciRef.names;

            switch bandInfer
                case 'ap'
                    ci.gains = cat(1, imec.sourceDatasets.apGain);
                    rmin = arrayfun(@(src) src.apRange(1), imec.sourceDatasets);
                    rmax = arrayfun(@(src) src.apRange(2), imec.sourceDatasets);
                    ci.ranges = cat(2, rmin, rmax);
                    
                case 'lf'
                    ci.gains = cat(1, imec.sourceDatasets.lfGain);
                    rmin = arrayfun(@(src) src.lfRange(1), imec.sourceDatasets);
                    rmax = arrayfun(@(src) src.lfRange(2), imec.sourceDatasets);
                    ci.ranges = cat(2, rmin, rmax);
            end
            
            [ci.multipliers, ~] = Neuropixel.ImecDataset.determineCommonGain(ci.gains, true); % quiet here since the gain is inferred anyway
            ci.adcBits = cat(1, imec.sourceDatasets.adcBits);
        end
    end
    
    methods
        function n = get.nDatasets(ds)
            n = numel(ds.samplesPreShift);
        end
        
        function v = get.samples(ci)
            if isempty(ci.timeShifts)
                v = ci.samplesPreShift;
            else
                v = arrayfun(@(spec) spec.nShiftedTimes, ci.timeShifts);
            end
        end
        
        function v = get.startIdx(ci)
            cs = cumsum([uint64(1); ci.samples]);
            v = cs(1:end-1);
        end
        
        function v = get.stopIdx(ci)
            v = [ci.startIdx(2:end); ci.nSamples];
        end
        
        function v = get.scaleToUvs(ci)
            v = (ci.ranges(:, 2) - ci.ranges(:, 1)) ./ (2.^ci.adcBits) ./ ci.gains * 1e6;
        end
        
        function [fileInds, origSampleInds] = lookup_sampleIndexInSourceFiles(ci, inds)
           % convert from a time index in this file to which of the concatenated files that index came from
           % factoring in the time shifts that were used to excise regions from those individual files before concatenating

           [fileInds, origSampleInds] = deal(zeros(size(inds), 'uint64'));
           starts = ci.startIdx;
           stops = ci.stopIdx;

           for i = 1:numel(starts)
              mask = inds >= starts(i) & inds <= stops(i);
              fileInds(mask) = i;
              origSampleInds(mask) = inds(mask) - double(starts(i)) + 1;
           end
           
           % now unshift the original sampleInds
           if ~isempty(ci.timeShifts)
               for i = 1:ci.nDatasets
                   mask = fileInds == i;
                   origSampleInds(mask) = ci.timeShifts(i).unshiftTimes(origSampleInds(mask));
               end
           end
        end
        
        function markConcatenatedFileBoundaries(ci, varargin)
            % mark the time points where multiple files were concatenated together
            % this assumes we are plotting into a space where ci.timeShifts is already applied (i.e. the concatenated file time axis),
            % but time_shifts may be passed in if needed to simulate a second set of time shifts
            p = inputParser();
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % masks which trials to include, note that it applies before the time shifts
            p.addParameter('side', 'bottom', @ischar);
            p.addParameter('Color', [0.3 0.3 0.9], @isvector);
            p.addParameter('LineWidth', 1, @isscalar);
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            p.addParameter('expand_limits', false, @islogical);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            timeShifts = p.Results.time_shifts; %#ok<*PROPLC>
            
            starts = ci.startIdx;
            if ~isempty(timeShifts)
                starts = timeShifts.shiftTimes(starts);
            end
            if ~isempty(p.Results.sample_window)
                window = p.Results.sample_window;
                mask = starts >= window(1) & starts <= window(2);
            else
                mask = true(size(starts));
            end
            if p.Results.timeInSeconds
                starts = double(starts) / ci.fs;
            end
            washolding = ishold();
            for iF = 1:numel(starts)
                if ~mask(iF), continue, end
                h = xline(starts(iF) + xOffset, '-', ci.names{iF}, 'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth, 'Interpreter', 'none');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
                hold on;
            end
            if ~washolding, hold off; end
        end
        
        function markExcisionBoundaries(ci, varargin)
            % marks times where each concatenated source file was split based on ci.timeShifts(iF)
            % this assumes we are plotting into a space where ci.timeShifts is NOT already applied (i.e. the concatenated file time axis)
            % but if time_shifts is passed in as ci.time_shifts, it will plot into the already time-shifted time axis
            p = inputParser();
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x));
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % assumed to be applied to times in plot
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            sample_window = p.Results.sample_window;
            
            timeInSeconds = p.Results.timeInSeconds;
            if isempty(ci.timeShifts)
                return;
            end
            
            % map the sample_window into individual files
            if isempty(sample_window)
                sample_window = [1 ci.nSamples];
            end
            [window_fileInd, window_origInd] = ci.lookup_sampleIndexInSourceFiles(sample_window);
            
            for iF = 1:ci.nDatasets
                if window_fileInd(1) > iF || window_fileInd(2) < iF, continue, end % this file isn't in window at all
                
                sample_window_this = [-Inf Inf];
                if window_fileInd(1) == iF
                    sample_window_this(1) = window_origInd(1);
                end
                if window_fileInd(2) == iF
                    sample_window_this(2) = window_origInd(2);
                end
                ci.timeShifts(iF).markExcisionBoundaries('sample_window', sample_window_this, 'xOffset', xOffset, ...
                    'timeInSeconds', timeInSeconds, 'fs', ci.fs, 'time_shifts', p.Results.time_shifts);
                if timeInSeconds
                    xOffset = xOffset + ci.timeShifts(iF).nShiftedTimes / ci.fs;
                else
                    xOffset = xOffset + ci.timeShifts(iF).nShiftedTimes;
                end
            end
        end
    end
end
        