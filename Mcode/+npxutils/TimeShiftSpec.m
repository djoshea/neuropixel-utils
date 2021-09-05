classdef TimeShiftSpec < handle & matlab.mixin.Copyable
    % Specification of a time shift.
    %
    % A simple class for holding information about a time shift that takes a set
    % of specific time intervals in the original file and maps them to new
    % indices that have been squished together by excising some intervening
    % intervals This is used primarily for skipping over bad regions in a raw
    % data file we don't want to sort
    
    properties
        % Original start index for interval.
        idxStart (:,1) uint64
        % Last index within interval.
        idxStop (:,1) uint64
        % New, shifted start index for interval.
        idxShiftStart (:,1) uint64
    end
    
    properties (Dependent)
        % Number of intervals.
        nIntervals
        % Durations of each interval.
        intervalDurations
        % Number of output samples that the source will be compressed into.
        nShiftedTimes
    end
    
    methods
        function this = TimeShiftSpec(idxStart, idxStop, idxShiftStart)
            % Construct a new object
            %
            % npxutils.TimeShiftSpec(idxStart, idxStop, idxShiftStart)
            % npxutils.TimeShiftSpec(timeShiftSpec)
            if nargin == 1
                if isa(idxStart, 'npxutils.TimeShiftSpec')
                    this = idxStart;
                else
                    mat = idxStart;
                    this.idxStart = uint64(mat(:, 1));
                    this.idxStop = uint64(mat(:, 2));
                    this.idxShiftStart = uint64(mat(:, 3));
                end
            elseif nargin == 3
                this.idxStart = uint64(idxStart);
                this.idxStop = uint64(idxStop);
                this.idxShiftStart = uint64(idxShiftStart);
            end
        end
        
        function specNew = convertToDifferentSampleRate(this, fsThis, fsNew, maxSamples)
            if nargin < 4
                maxSamples = Inf;
            end
            clamp = @(idx) min(max(idx, 1), maxSamples);
            convert = @(idx) uint64(clamp(floor(double(idx) / double(fsThis) * double(fsNew))));
            specNew = npxutils.TimeShiftSpec(convert(this.idxStart), convert(this.idxStop), convert(this.idxShiftStart));
        end
        
        function mat = as_matrix(this)
            % Convert to a matrix.
            %
            % See also:
            % from_matrixd
            mat = cat(2, this.idxStart, this.idxStop, this.idxShiftStart);
        end
        
        function str = as_string(this)
            % Convert to a string.
            %
            % See also:
            % from_string
            mat = this.as_matrix();
            if isempty(mat)
                str = "[]";
            else
                str = string(mat2str(mat));
            end
        end
        
        function dur = get.intervalDurations(this)
            dur = uint64(int64(this.idxStop) - int64(this.idxStart)) + uint64(1);
            dur(this.idxStart == uint64(0) | this.idxStop == uint64(0)) = 0;
        end
        
        function n = get.nIntervals(this)
            n = numel(this.idxStart);
        end
        
        function n = get.nShiftedTimes(this)
            lastOutSample = uint64(int64(this.idxShiftStart) + int64(this.intervalDurations) - int64(1));
            n = max(lastOutSample);
        end
        
        function [times_out, mask] = shiftTimes(this, times)
            % Reassign given times by shifting them.
            %
            % [times_out, mask] = shiftTimes(this, times)
            %
            % Reassigns the times provided by shifting them according the delta
            % between idxStart and idxShiftStart. If a time doesn't fit within
            % any of the specified windows, it will be marked as NaN.
            %
            % Returns:
            %   times_out - the shifted times.
            %   mask - a logical vector indicating whether each time fit within
            %     a window.
            %
            % See also:
            % unshiftTimes
            mask = false(numel(times), 1);
            times_out = times;
            for iR = 1:this.nIntervals
                match_this = times >= this.idxStart(iR) & times <= this.idxStop(iR);
                mask(match_this) = true;
                times_out(match_this) = uint64(int64(times(match_this)) - int64(this.idxStart(iR)) + int64(this.idxShiftStart(iR)));
            end
            
            times_out(~mask) = NaN;
        end
        
        function time_src = unshiftTimes(this, times)
            % Recover the original time that was mapped to times by shiftTimes
            %
            % time_src = unshiftTimes(this, times)
            %
            % Times is the shifted times you want to unshift.
            %
            % See also:
            % shiftTimes
            dur = this.intervalDurations;
            
            time_src = zeros(size(times), 'like', times);
            for iR = 1:this.nIntervals
                mask = times >= this.idxShiftStart(iR) & times <= this.idxShiftStart(iR) + dur(iR); % times in this interval
                % no offset by 1 here, time == idxShiftStart should map to idxStart exactly
                time_src(mask) = uint64(int64(times(mask)) - int64(this.idxShiftStart(iR)) +  int64(this.idxStart(iR)));
            end
            %             srcIdx = spec.computeSourceIndices();
            %             [~, time_src] = ismember(times, srcIdx);
        end
        
        function constrainToNumSamples(this, nSamplesSource)
            maskKeep = this.idxStart <= uint64(nSamplesSource);
            this.idxStart = this.idxStart(maskKeep);
            this.idxStop = min(this.idxStop(maskKeep), uint64(nSamplesSource));
            this.idxShiftStart = this.idxShiftStart(maskKeep);
        end
        
        function constrainToSampleWindow(this, window)
            maskKeep = this.idxStop >= uint64(window(1)) & this.idxStart <= uint64(window(2));
            this.idxStart = max(uint64(window(1)), this.idxStart(maskKeep));
            this.idxStop = min(this.idxStop(maskKeep), uint64(window(2)));
            this.idxShiftStart = this.idxShiftStart(maskKeep);
        end
        
        function shiftIndices = computeSourceIndices(this, nSamplesSource)
            % Compute source IMEC indices for given number of samples.
            %
            % Computes a nShiftedTimes x 1 vector of indices into the source
            % file imec that would occupy times 1:nShiftedTimes in the shifted
            % output file. If a given slot isn't filled, it will be set to 0.
            
            if nargin > 1
                this.constrainToNumSamples(nSamplesSource);
            end
            
            shiftIndices = zeros(this.nShiftedTimes, 1, 'uint64');
            dur = this.intervalDurations;
            for i = 1:this.nIntervals
                offsets = int64(0):(int64(dur(i)) - int64(1));
                to = int64(this.idxShiftStart(i)) + offsets;
                from = int64(this.idxStart(i)) + offsets;
                shiftIndices(to) = from;
            end
        end
        
        function h = markExcisionBoundaries(this, varargin)
            p = inputParser();
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % set this to the x window being shown to avoid plotting hidden boundaries
            p.addParameter('Color', [0.9 0.3 0.9], @isvector);
            p.addParameter('LineWidth', 1, @isscalar);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples
            p.addParameter('fs', 0, @(x) isempty(x) || isscalar(x));
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec')); % applied on what is plotted, on top of this spec
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            boundaries = double(this.idxShiftStart(2:end));
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
    
    methods (Static)
        function spec = non_shift(nSamples)
            % Create a TimeShiftSpec that takes all samples without any shifting.
            spec = npxutils.TimeShiftSpec(1, nSamples, 1);
        end
        
        function spec = take_no_samples()
            % Create a TimeShiftSpec that takes no samples.
            spec = npxutils.TimeShiftSpec(0, 0, 1);
        end
        
        function spec = from_matrix(mat, nSamplesFull)
            % Create a TimeShiftSpec from a matrix.
            %
            % spec = from_matrix(mat, nSamplesFull)
            %
            % mat is cat(2, spec.idxStart, spec.idxStop, spec.idxShiftStart).
            if isempty(mat)
                spec = npxutils.TimeShiftSpec.non_shift(nSamplesFull);
            else
                assert(size(mat, 2) == 3);
                spec = npxutils.TimeShiftSpec(mat(:, 1), mat(:, 2), mat(:, 3));
            end
        end
        
        function spec = from_string(str, nSamplesFull)
            % Create a TimeShiftSpec from a string representation.
            %
            % spec = from_string(str, nSamplesFull)
            %
            % Str is as generated from to_string, possibly multiple strings
            % concatenated.
            str = strtrim(string(str));
            if str == ""
                spec = npxutils.TimeShiftSpec.non_shift(nSamplesFull);
            elseif str == "[]"
                spec = npxutils.TimeShiftSpec.take_no_samples();
            else
                mat = str2num(str); %#ok<ST2NM>
                spec = npxutils.TimeShiftSpec.from_matrix(mat, nSamplesFull);
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
            
            spec = npxutils.TimeShiftSpec(idxStart, idxStop, regionUpdatedStarts);
        end
        
        function spec = buildToConstrainSampleWindow(window)
            spec = npxutils.TimeShiftSpec(window(1), window(2), uint64(1));
        end
        
        function spec = buildToVisualizeMultipleWindows(idxStart, idxStop, varargin)
            % Create a TimeShiftSpec to visualize multiple windows.
            %
            % spec = buildToVisualizeMultipleWindows(idxStart, idxStop, ...)
            %
            % Options:
            %   padding - defaults to 100.
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
            
            spec = npxutils.TimeShiftSpec(idxStart, idxStop, regionUpdatedStarts);
        end
    end
    
end