classdef TrialSegmentationInfo < handle & matlab.mixin.Copyable
    % simple utility for describing how to segment trials in Kilosort

    properties
        fs % sampling rate (Hz)
        trialId(:, 1) uint32
        conditionId(:, 1) uint32
        idxStart(:, 1) uint64 % in samples
        idxStop(:, 1) uint64 % in samples, may not be used
    end
    
    properties(Dependent)
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
        
        function tsiNew = convertToDifferentSampleRate(tsi, fsNew)
            tsiNew = copy(tsi);
            tsiNew.fs = fsNew;
            
            convert = @(idx) uint64(floor(double(idx) / double(tsi.fs) * double(fsNew)));
            tsiNew.idxStart = convert(tsi.idxStart);
            tsiNew.idxStop = convert(tsi.idxStop);
        end
        
        function nTrials = get.nTrials(tsi)
            nTrials = numel(tsi.trialId);
        end
        
        function maskTrials(tsi, mask)
            tsi.trialId = tsi.trialId(mask);
            tsi.conditionId = tsi.conditionId(mask);
            tsi.idxStart = tsi.idxStart(mask);
            tsi.idxStop = tsi.idxStop(mask);
        end
        
        function maskTrialsByTrialId(tsi, trialIds)
            mask = ismember(tsi.trialId, trialIds);
            tsi.maskTrials(mask);
        end
        
        function maskTrialsWithinTrialIdRange(tsi, trialIdLims)
            trialIdLims = uint32(trialIdLims);
            assert(numel(trialIdLims) == 2);
            mask = tsi.trialId >= trialIdLims(1) & tsi.trialId <= trialIdLims(2);
            tsi.maskTrials(mask);
        end

        function maskTimeWindowMsWithinTrial(tsi, offsetMsStart, offsetMsStop)
            % reduces the size of each trial to a specific time window within each trial relative to start
            assert(numel(offsetMsStart) == tsi.nTrials);
            assert(numel(offsetMsStop) == tsi.nTrials);

            offsetSampleStart = int64(round(double(offsetMsStart) * double(tsi.fs) / 1000));
            offsetSampleStop = int64(round(double(offsetMsStop) * double(tsi.fs) / 1000));
            start = int64(tsi.idxStart);
            tsi.idxStart = uint64(start + offsetSampleStart);
            tsi.idxStop = uint64(start + offsetSampleStop);
        end
        
        function tf = get.trialsAreAdjacent(tsi)
            % true if each trial's stop == next trial's start (or within one sample
            tf = all(int64(tsi.idxStart(2:end)) - int64(tsi.idxStop(1:end-1)) <= int64(1));
        end
        
        function [trialInds, trialIds, mask_in_trial] = segmentTimes(tsi, times)
            edges = [tsi.idxStart; tsi.idxStop(end) + uint64(1)];
            trialInds = discretize(times, edges);
            trialInds(trialInds == numel(edges)) = NaN; % ignore the last bin
            if ~tsi.trialsAreAdjacent
                % nan out trial assignments for times which lie past the corresponding trial stop
                mask = ~isnan(trialInds);
                correspondingStop = zeros(size(times), 'uint64');
                correspondingStop(mask) = tsi.idxStop(trialInds(mask));
                trialInds(times > correspondingStop) = NaN;
            end
            
            trialIds = zeros(size(times), 'like', trialInds);
            mask_in_trial = ~isnan(trialInds);
            trialIds(mask_in_trial) = tsi.trialId(trialInds(mask_in_trial));
        end
        
        function truncateTrialsLongerThan(tsi, maxDurSec)
            maxDurSamples = maxDurSec * tsi.fs;
            durations = tsi.idxStop - tsi.idxStart;
            longTrials = durations > maxDurSamples;
            tsi.idxStop(longTrials) = tsi.idxStart(longTrials) + maxDurSamples;
        end
     
        function [idxStart, idxStop, trialStartStop] = computeActiveRegions(tsi, varargin)
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar); % in samples
            p.parse(varargin{:});
            
            if tsi.nTrials == 0
                [idxStart, idxStop, trialStartStop] = deal(zeros(0, 1, 'uint64'));
                return;
            end
            
            % find regions where there are no pauses greater than maxPauseSec
            maxPauseSamples = p.Results.maxPauseSec * tsi.fs;
            pauses = int64(tsi.idxStart(2:end)) - int64(tsi.idxStop(1:end-1));
            longPauses = find(pauses > maxPauseSamples); % longPause after trial idx
            
            R = numel(longPauses);
            [idxStart, idxStop] = deal(zeros(R+1, 1, 'uint64'));
            trialStartStop = zeros(R+1, 2, 'uint32');
            last = 1;
            for iR = 1:R
                idxStart(iR) = tsi.idxStart(last);
                idxStop(iR) = tsi.idxStop(longPauses(iR));
                trialStartStop(iR, :) = [tsi.trialId(last), tsi.trialId(longPauses(iR))];
                last = longPauses(iR) + 1;
            end

            idxStart(end) = tsi.idxStart(last)-uint64(1);
            idxStop(end) = tsi.idxStop(end);
            trialStartStop(end, :) = [tsi.trialId(last), tsi.trialId(end)]; 
        end
             
        function markTrialStartStops(tsi, varargin) % assumes x is in seconds
            p = inputParser();
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            hold on;
            for iT = 1:numel(tsi.idxStart)
                xline(tsi.idxStart(iT) / tsi.fs + xOffset, p.Results.LineSpecStart, num2str(tsi.trialId(iT)), 'Color', p.Results.Color, 'LineWidth', 0.5);
                xline(tsi.idxStart(iT) / tsi.fs + xOffset, p.Results.LineSpecStop, 'Color', p.Results.Color, 'LineWidth', 0.5);
            end
        end
        
        function h = plotTrialsActiveRegions(tsi, varargin)
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar);
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            cla;
            ylim([0 1]);
            tsi.markActiveRegions(p.Results, p.Unmatched);
            tsi.markTrialTicks(p.Unmatched);
            
        end
        
        function h = markTrialTicks(tsi, varargin)
            p = inputParser();
            p.addParameter('maskTrials', true(tsi.nTrials, 1), @isvector);
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % masks which trials to include, note that it applies before the time shifts
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('side', 'bottom', @ischar);
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            p.addParameter('expand_limits', false, @islogical);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            timeShifts = p.Results.time_shifts;
            ticks =  tsi.idxStart;
            
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
                ticks = double(ticks) / tsi.fs;
            end
            
            templateRows = [ ...
                dataTipTextRow('idx start', tsi.idxStart(mask), '%d'); ...
                dataTipTextRow('trial id', tsi.trialId(mask), '%d'); ...
                dataTipTextRow('idx start', tsi.conditionId(mask), '%d') ];
            
            h = Neuropixel.Utils.rugplot(ticks + xOffset, 'side', p.Results.side, 'Color', p.Results.Color, 'expand_limits', p.Results.expand_limits, ...
                'dataTipTemplateRows', templateRows);
            set(h, 'XLimInclude', 'off');
        end
        
        function markActiveRegions(tsi, varargin)
            % marks regions that are densely covered with trials (ignoring gaps between trials
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar); 
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            [idxStart, idxStop, trialStartStop] = tsi.computeActiveRegions('maxPauseSec', p.Results.maxPauseSec); %#ok<*PROPLC>
            timeShiftSpec = p.Results.time_shifts;
            if ~isempty(timeShiftSpec)
                idxStart = timeShiftSpec.shiftTimes(idxStart);
                idxStop = timeShiftSpec.shiftTimes(idxStop);
            end
            
            if p.Results.timeInSeconds
                idxStart = double(idxStart) / tsi.fs;
                idxStop = double(idxStop) / tsi.fs;
            end
            hold on;
            for iT = 1:numel(idxStart)
                str = sprintf('Trial %d-%d', trialStartStop(iT, 1), trialStartStop(iT, 2));
                xline(double(idxStart(iT)) + xOffset, p.Results.LineSpecStart, '', 'Color', p.Results.Color, 'LineWidth', 0.5, 'Interpreter', 'none');
                h = xline(double(idxStop(iT)) + xOffset, p.Results.LineSpecStop, str, 'Color', p.Results.Color, 'LineWidth', 0.5, 'Interpreter', 'none');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
            end
        end
        
        function markSpecificTrials(tsi, trialIds, varargin)
            % marks regions that are densely covered with trials (ignoring gaps between trials
            p = inputParser();
            p.addParameter('labels', [], @(x) isempty(x) || isstringlike(x));
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % masks which trials to include, note that it applies before the time shifts
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('Color', [0 0 1], @ismatrix);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            [tf, trialInds] = ismember(trialIds, tsi.trialId);
            assert(all(tf), "Some trialIds not found in TrialSegmentationInfo");
            
            labels = p.Results.labels;
            if isempty(labels)
                labels = arrayfun(@(t) sprintf("Trial %d", t), trialIds);
            else
                labels = string(labels);
            end
            
            idxStart = tsi.idxStart(trialInds);
            
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
                idxStart = double(idxStart) / tsi.fs;
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
        
        function spec = computeShiftsExciseRegionsOutsideTrials(tsi, varargin)
            % shifts is an nRegions x 3 matrix of uint64 sample idxs
            % start of original window, stop of original window, updated start for this window
            [idxStart, idxStop, ~] = tsi.computeActiveRegions(varargin{:}); %#ok<*PROP>
            spec = Neuropixel.TimeShiftSpec.buildToExciseGaps(idxStart, idxStop);
        end
    end

end
