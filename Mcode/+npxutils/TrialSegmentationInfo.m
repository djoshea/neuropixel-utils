classdef TrialSegmentationInfo < handle & matlab.mixin.Copyable
    % Simple utility for describing how to segment trials in Kilosort.
    
    properties
        % sampling rate in Hz
        fs
        trialId (:,1) uint32
        conditionId (:,1) uint32
        idxStart (:,1) uint64 % in samples
        idxStop (:,1) uint64 % in samples, may not be used
    end
    
    properties (Dependent)
        nTrials
        trialsAreAdjacent % true if each trial's stop == next trial's start
    end
    
    methods
        function tsi = TrialSegmentationInfo(nTrials, fs)
            if nargin < 1
                nTrials = 0;
            end
            if nargin < 2
                fs = [];
            end
            tsi.fs = fs;
            tsi.trialId = zeros(nTrials, 1, 'uint32');
            tsi.conditionId = zeros(nTrials, 1, 'uint32');
            tsi.idxStart = zeros(nTrials, 1, 'uint64');
            tsi.idxStop = zeros(nTrials, 1, 'uint64');
        end
        
        function tsiNew = convertToDifferentSampleRate(this, fsNew)
            tsiNew = copy(this);
            tsiNew.fs = fsNew;
            
            convert = @(idx) uint64(floor(double(idx) / double(this.fs) * double(fsNew)));
            tsiNew.idxStart = convert(this.idxStart);
            tsiNew.idxStop = convert(this.idxStop);
        end
        
        function nTrials = get.nTrials(this)
            nTrials = numel(this.trialId);
        end
        
        function maskTrials(this, mask)
            this.trialId = this.trialId(mask);
            this.conditionId = this.conditionId(mask);
            this.idxStart = this.idxStart(mask);
            this.idxStop = this.idxStop(mask);
        end
        
        function maskTrialsByTrialId(this, trialIds)
            mask = ismember(this.trialId, trialIds);
            this.maskTrials(mask);
        end
        
        function maskTrialsWithinTrialIdRange(this, trialIdLims)
            trialIdLims = uint32(trialIdLims);
            assert(numel(trialIdLims) == 2);
            mask = this.trialId >= trialIdLims(1) & this.trialId <= trialIdLims(2);
            this.maskTrials(mask);
        end
        
        function tf = get.trialsAreAdjacent(this)
            % true if each trial's stop == next trial's start (or within one sample
            tf = all(int64(this.idxStart(2:end)) - int64(this.idxStop(1:end-1)) <= int64(1));
        end
        
        function [trialInds, trialIds, mask_in_trial] = segmentTimes(this, times)
            edges = [this.idxStart; this.idxStop(end) + uint64(1)];
            trialInds = discretize(times, edges);
            trialInds(trialInds == numel(edges)) = NaN; % ignore the last bin
            if ~this.trialsAreAdjacent
                % nan out trial assignments for times which lie past the corresponding trial stop
                mask = ~isnan(trialInds);
                correspondingStop = zeros(size(times), 'uint64');
                correspondingStop(mask) = this.idxStop(trialInds(mask));
                trialInds(times > correspondingStop) = NaN;
            end
            
            trialIds = zeros(size(times), 'like', trialInds);
            mask_in_trial = ~isnan(trialInds);
            trialIds(mask_in_trial) = this.trialId(trialInds(mask_in_trial));
        end
        
        function truncateTrialsLongerThan(this, maxDurSec)
            maxDurSamples = maxDurSec * this.fs;
            durations = this.idxStop - this.idxStart;
            longTrials = durations > maxDurSamples;
            this.idxStop(longTrials) = this.idxStart(longTrials) + maxDurSamples;
        end
        
        function [idxStart, idxStop, trialStartStop] = computeActiveRegions(this, varargin)
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar); % in samples
            p.parse(varargin{:});
            
            if this.nTrials == 0
                [idxStart, idxStop, trialStartStop] = deal(zeros(0, 1, 'uint64'));
                return;
            end
            
            % find regions where there are no pauses greater than maxPauseSec
            maxPauseSamples = p.Results.maxPauseSec * this.fs;
            pauses = int64(this.idxStart(2:end)) - int64(this.idxStop(1:end-1));
            longPauses = find(pauses > maxPauseSamples); % longPause after trial idx
            
            R = numel(longPauses);
            [idxStart, idxStop] = deal(zeros(R+1, 1, 'uint64'));
            trialStartStop = zeros(R+1, 2, 'uint32');
            last = 1;
            for iR = 1:R
                idxStart(iR) = this.idxStart(last);
                idxStop(iR) = this.idxStop(longPauses(iR));
                trialStartStop(iR, :) = [this.trialId(last), this.trialId(longPauses(iR))];
                last = longPauses(iR) + 1;
            end
            
            idxStart(end) = this.idxStart(last)-uint64(1);
            idxStop(end) = this.idxStop(end);
            trialStartStop(end, :) = [this.trialId(last), this.trialId(end)];
        end
        
        function markTrialStartStops(this, varargin) % assumes x is in seconds
            p = inputParser();
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            hold on;
            for iT = 1:numel(this.idxStart)
                xline(this.idxStart(iT) / this.fs + xOffset, p.Results.LineSpecStart, num2str(this.trialId(iT)), 'Color', p.Results.Color, 'LineWidth', 0.5);
                xline(this.idxStart(iT) / this.fs + xOffset, p.Results.LineSpecStop, 'Color', p.Results.Color, 'LineWidth', 0.5);
            end
        end
        
        function h = plotTrialsActiveRegions(this, varargin)
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar);
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            cla;
            ylim([0 1]);
            this.markActiveRegions(p.Results, p.Unmatched);
            this.markTrialTicks(p.Unmatched);
            
        end
        
        function h = markTrialTicks(this, varargin)
            p = inputParser();
            p.addParameter('maskTrials', true(this.nTrials, 1), @isvector);
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % masks which trials to include, note that it applies before the time shifts
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec'));
            p.addParameter('side', 'bottom', @ischar);
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples
            p.addParameter('expand_limits', false, @islogical);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            timeShifts = p.Results.time_shifts;
            ticks =  this.idxStart;
            
            if ~isempty(timeShifts)
                ticks = timeShifts.shiftTimes(ticks);
            end
            sample_window = p.Results.sample_window;
            if ~isempty(sample_window)
                mask = ticks >= sample_window(1) & ticks <= sample_window(2);
            else
                mask = true(size(ticks));
            end
            mask = mask & p.Results.maskTrials;
            ticks = ticks(mask);
            
            if p.Results.timeInSeconds
                ticks = double(ticks) / this.fs;
            end
            
            templateRows = [ ...
                dataTipTextRow('idx start', this.idxStart(mask), '%d'); ...
                dataTipTextRow('trial id', this.trialId(mask), '%d'); ...
                dataTipTextRow('idx start', this.conditionId(mask), '%d') ];
            
            h = npxutils.internal.graphics.rugplot(ticks + xOffset, 'side', p.Results.side, 'Color', p.Results.Color, 'expand_limits', p.Results.expand_limits, ...
                'dataTipTemplateRows', templateRows);
            set(h, 'XLimInclude', 'off');
        end
        
        function markActiveRegions(this, varargin)
            % marks regions that are densely covered with trials (ignoring gaps between trials
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar);
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec'));
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples
            
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            [idxStart, idxStop, trialStartStop] = this.computeActiveRegions('maxPauseSec', p.Results.maxPauseSec); %#ok<*PROPLC>
            timeShiftSpec = p.Results.time_shifts;
            if ~isempty(timeShiftSpec)
                idxStart = timeShiftSpec.shiftTimes(idxStart);
                idxStop = timeShiftSpec.shiftTimes(idxStop);
            end
            
            if p.Results.timeInSeconds
                idxStart = double(idxStart) / this.fs;
                idxStop = double(idxStop) / this.fs;
            end
            hold on;
            for iT = 1:numel(idxStart)
                str = sprintf('Trial %d-%d', trialStartStop(iT, 1), trialStartStop(iT, 2));
                xline(double(idxStart(iT)) + xOffset, p.Results.LineSpecStart, '', 'Color', p.Results.Color, 'LineWidth', 0.5, 'Interpreter', 'none');
                h = xline(double(idxStop(iT)) + xOffset, p.Results.LineSpecStop, str, 'Color', p.Results.Color, 'LineWidth', 0.5, 'Interpreter', 'none');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
            end
        end
        
        function markSpecificTrials(this, trialIds, varargin)
            % marks regions that are densely covered with trials (ignoring gaps between trials
            p = inputParser();
            p.addParameter('labels', [], @(x) isempty(x) || isstringlike(x));
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % masks which trials to include, note that it applies before the time shifts
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec'));
            p.addParameter('Color', [0 0 1], @ismatrix);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples
            
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            [tf, trialInds] = ismember(trialIds, this.trialId);
            assert(all(tf), "Some trialIds not found in TrialSegmentationInfo");
            
            labels = p.Results.labels;
            if isempty(labels)
                labels = arrayfun(@(t) sprintf("Trial %d", t), trialIds);
            else
                labels = string(labels);
            end
            
            idxStart = this.idxStart(trialInds);
            
            sample_window = p.Results.sample_window;
            if ~isempty(sample_window)
                mask = idxStart >= sample_window(1) & idxStart <= sample_window(2);
            else
                mask = true(size(idxStart));
            end
            idxStart = idxStart(mask);
            
            timeShiftSpec = p.Results.time_shifts;
            if ~isempty(timeShiftSpec)
                idxStart = timeShiftSpec.shiftTimes(idxStart);
            end
            
            if p.Results.timeInSeconds
                idxStart = double(idxStart) / this.fs;
            end
            hold on;
            
            cmap = p.Results.Color;
            if size(cmap, 1) == 1
                cmap = repmat(cmap, numel(idxStart), 1);
            else
                cmap = cmap(mask, :);
            end
            for iT = 1:numel(idxStart)
                h = xline(double(idxStart(iT)) + xOffset, p.Results.LineSpecStart, labels{iT}, 'Color', cmap(iT, :), 'LineWidth', 0.5, 'Interpreter', 'none');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
            end
        end
        
        function spec = computeShiftsExciseRegionsOutsideTrials(this, varargin)
            % shifts is an nRegions x 3 matrix of uint64 sample idxs
            % start of original window, stop of original window, updated start for this window
            [idxStart, idxStop, ~] = this.computeActiveRegions(varargin{:}); %#ok<*PROP>
            spec = npxutils.TimeShiftSpec.buildToExciseGaps(idxStart, idxStop);
        end
    end
    
end
