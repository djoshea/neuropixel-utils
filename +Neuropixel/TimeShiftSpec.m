classdef TimeShiftSpec < handle & matlab.mixin.Copyable
    % a simple class for holding information about a time shift that takes a set of specific time intervals
    % in the original file and maps them to new indices that have been squished together by excising some intervening intervals
    % This is used primarily for skipping over bad regions in a raw data file we don't want to sort
    
    properties
        idxStart(:, 1) uint64 % original start index for interval
        idxStop(:, 1) uint64 % last index within interval
        idxShiftStart(:, 1) uint64 % new, shifted start index for interval
    end
    
    properties(Dependent)
        nIntervals
        intervalDurations
        nShiftedTimes % number of output samples that the source will be compressed into
    end
    
    methods
        function spec = TimeShiftSpec(idxStart, idxStop, idxShiftStart)
            if nargin == 1
                if isa(idxStart, 'Neuropixel.TimeShiftSpec')
                    spec = idxStart;
                else
                    mat = idxStart;
                    spec.idxStart = uint64(mat(:, 1));
                    spec.idxStop = uint64(mat(:, 2));
                    spec.idxShiftStart = uint64(mat(:, 3));
                end
            elseif nargin == 3
                spec.idxStart = uint64(idxStart);
                spec.idxStop = uint64(idxStop);
                spec.idxShiftStart = uint64(idxShiftStart);
            end     
        end
        
        function specNew = convertToDifferentSampleRate(spec, fsThis, fsNew, maxSamples)
            if nargin < 4
                maxSamples = Inf;
            end
            clamp = @(idx) min(max(idx, 1), maxSamples);
            convert = @(idx) uint64(clamp(floor(double(idx) / double(fsThis) * double(fsNew))));
            specNew = Neuropixel.TimeShiftSpec(convert(spec.idxStart), convert(spec.idxStop), convert(spec.idxShiftStart));
        end
        
        function mat = as_matrix(spec)
            % see from_matrixd
            mat = cat(2, spec.idxStart, spec.idxStop, spec.idxShiftStart);
        end
        
        function str = as_string(spec)
            % see from_string
            mat = spec.as_matrix();
            if isempty(mat)
                str = "[]";
            else
                str = string(mat2str(mat));
            end
        end
        
        function dur = get.intervalDurations(spec)
            dur = uint64(int64(spec.idxStop) - int64(spec.idxStart)) + uint64(1);
            dur(spec.idxStart == uint64(0) | spec.idxStop == uint64(0)) = 0;
        end
        
        function n = get.nIntervals(spec)
            n = numel(spec.idxStart);
        end
        
        function n = get.nShiftedTimes(spec)
            lastOutSample = uint64(int64(spec.idxShiftStart) + int64(spec.intervalDurations) - int64(1));
            n = max(lastOutSample);
        end
        
        function [times_out, mask] = shiftTimes(spec, times)
            % reassign the times provided by shifting them according the delta between idxStart and idxShiftStart
            % if a time doesn't fit within any of the specified windows, it will be marked as NaN.
            mask = false(numel(times), 1);
            times_out = times;
            for iR = 1:spec.nIntervals
                match_this = times >= spec.idxStart(iR) & times <= spec.idxStop(iR);
                mask(match_this) = true;
                times_out(match_this) = uint64(int64(times(match_this)) - int64(spec.idxStart(iR)) + int64(spec.idxShiftStart(iR)));
            end
            
            times_out(~mask) = NaN;
        end
        
        function time_src = unshiftTimes(spec, times)
            % recover the original time that was mapped to times by shiftTimes
            dur = spec.intervalDurations;
            
            time_src = zeros(size(times), 'like', times);
            for iR = 1:spec.nIntervals
                mask = times >= spec.idxShiftStart(iR) & times <= spec.idxShiftStart(iR) + dur(iR); % times in this interval
                % no offset by 1 here, time == idxShiftStart should map to idxStart exactly
                time_src(mask) = uint64(int64(times(mask)) - int64(spec.idxShiftStart(iR)) +  int64(spec.idxStart(iR)));
            end
%             srcIdx = spec.computeSourceIndices();
%             [~, time_src] = ismember(times, srcIdx);
        end
        
        function constrainToNumSamples(spec, nSamplesSource)
            maskKeep = spec.idxStart <= uint64(nSamplesSource);
            spec.idxStart = spec.idxStart(maskKeep);
            spec.idxStop = min(spec.idxStop(maskKeep), uint64(nSamplesSource));
            spec.idxShiftStart = spec.idxShiftStart(maskKeep);
        end
        
        function constrainToSampleWindow(spec, window)
            maskKeep = spec.idxStop >= uint64(window(1)) & spec.idxStart <= uint64(window(2));
            spec.idxStart = max(uint64(window(1)), spec.idxStart(maskKeep));
            spec.idxStop = min(spec.idxStop(maskKeep), uint64(window(2)));
            spec.idxShiftStart = spec.idxShiftStart(maskKeep);
        end
        
        function shiftIndices = computeSourceIndices(spec, nSamplesSource)
            % compute a nShiftedTimes x 1 vector of indices into the source file imec
            % that would occupy times 1:nShiftedTimes in the shifted output file
            % if a given slot isn't filled, it will be set to 0
           
            if nargin > 1
                spec.constrainToNumSamples(nSamplesSource);
            end

            shiftIndices = zeros(spec.nShiftedTimes, 1, 'uint64');
            dur = spec.intervalDurations;
            for i = 1:spec.nIntervals
               offsets = int64(0):(int64(dur(i)) - int64(1));
               to = int64(spec.idxShiftStart(i)) + offsets;
               from = int64(spec.idxStart(i)) + offsets;
               shiftIndices(to) = from;
            end
        end
        
        function h = markExcisionBoundaries(spec, varargin)
            p = inputParser();
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % set this to the x window being shown to avoid plotting hidden boundaries
            p.addParameter('Color', [0.9 0.3 0.9], @isvector);
            p.addParameter('LineWidth', 1, @isscalar);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            p.addParameter('fs', 0, @(x) isempty(x) || isscalar(x));
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % applied on what is plotted, on top of this spec
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            boundaries = double(spec.idxShiftStart(2:end));
            time_shifts_plot = p.Results.time_shifts;
            if ~isempty(time_shifts_plot)
                boundaries = time_shifts_plot.shiftTimes(boundaries);
            end
            sample_window = p.Results.sample_window;
            if ~isempty(sample_window)
                mask = boundaries >= sample_window(1) & boundaries <= sample_window(2);
            else
                mask = true(size(boundaries));
            end
            if p.Results.timeInSeconds
                boundaries = boundaries / p.Results.fs;
            end
            washolding = ishold();
            for iF = 1:numel(boundaries)
                if ~mask(iF), continue, end
                h = xline(boundaries(iF) + xOffset, '-', sprintf('TimeShift %d', iF), 'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth, 'Interpreter', 'none');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
                hold on;
            end
            if ~washolding, hold off, end
        end
    end
    
    methods(Static)
        function spec = non_shift(nSamples)
            % returns a TimeShiftSpec that takes all samples without any shifting
            spec = Neuropixel.TimeShiftSpec(1, nSamples, 1);
        end
        
        function spec = take_no_samples()
            % returns a TimeShiftSpec that takes no samples
            spec = Neuropixel.TimeShiftSpec(0, 0, 1);
        end
        
        function spec = from_matrix(mat, nSamplesFull)
            % mat is mat = cat(2, spec.idxStart, spec.idxStop, spec.idxShiftStart)
            if isempty(mat)
                spec = Neuropixel.TimeShiftSpec.non_shift(nSamplesFull);
            else
                assert(size(mat, 2) == 3);
                spec = Neuropixel.TimeShiftSpec(mat(:, 1), mat(:, 2), mat(:, 3));
            end
        end
        
        function spec = from_string(str, nSamplesFull)
            % string is as generated from to_string, possibly multiple strings concatenated 
            str = strtrim(string(str));
            if str == ""
                spec = Neuropixel.TimeShiftSpec.non_shift(nSamplesFull);
            elseif str == "[]"
                spec = Neuropixel.TimeShiftSpec.take_no_samples();
            else
                mat = str2num(str); %#ok<ST2NM>
                spec = Neuropixel.TimeShiftSpec.from_matrix(mat, nSamplesFull);
            end
        end  
        
        function spec = buildToExciseGaps(idxStart, idxStop)
            regionDurations = int64(idxStop) - int64(idxStart);
            regionUpdatedStarts = int64(idxStart);
            if ~isempty(regionUpdatedStarts)
                regionUpdatedStarts(1) = 1;
                for iR = 2:numel(idxStart)
                    regionUpdatedStarts(iR) = regionUpdatedStarts(iR-1) + regionDurations(iR-1);
                end
            end

            spec = Neuropixel.TimeShiftSpec(idxStart, idxStop, regionUpdatedStarts);
        end
        
        function spec = buildToConstrainSampleWindow(window)
            spec = Neuropixel.TimeShiftSpec(window(1), window(2), uint64(1));
        end
        
        function spec = buildToVisualizeMultipleWindows(idxStart, idxStop, varargin)
            p = inputParser();
            p.addParameter('padding', 100, @isscalar);
            p.parse(varargin{:});
            
            padding = p.Results.padding;
            
            regionDurations = int64(idxStop) - int64(idxStart);
            regionUpdatedStarts = int64(idxStart);
            regionUpdatedStarts(1) = 1;
            for iR = 2:numel(idxStart)
                regionUpdatedStarts(iR) = regionUpdatedStarts(iR-1) + regionDurations(iR-1) + padding;
            end
            
            spec = Neuropixel.TimeShiftSpec(idxStart, idxStop, regionUpdatedStarts);
        end
    end
    
end