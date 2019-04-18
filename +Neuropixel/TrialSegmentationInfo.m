classdef TrialSegmentationInfo < handle & matlab.mixin.Copyable
    % simple utility for describing how to segment trials in KiloSort

    properties
        fs % sampling rate (Hz)
        trialId(:, 1) uint32
        conditionId(:, 1) uint32
        idxStart(:, 1) uint64 % in samples
        idxStop(:, 1) uint64 % in samples, may not be used
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
            p.parse(varargin{:});
            
            hold on;
            for iT = 1:numel(tsi.idxStart)
                xline(tsi.idxStart(iT) / tsi.fs, p.Results.LineSpecStart, num2str(tsi.trialId(iT)), 'Color', p.Results.Color, 'LineWidth', 0.5);
                xline(tsi.idxStart(iT) / tsi.fs, p.Results.LineSpecStop, 'Color', p.Results.Color, 'LineWidth', 0.5);
            end
        end
        
        function h = markTrialTicks(tsi, varargin)
            p = inputParser();
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % masks which trials to include, note that it applies before the time shifts
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('side', 'bottom', @ischar);
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('timeInSeconds', true, @islogical); % if true, x axis is seconds, if false, is samples 
            p.addParameter('expand_limits', false, @islogical);
            p.parse(varargin{:});
            
            timeShifts = p.Results.time_shifts;
            ticks =  tsi.idxStart;
            sample_window = p.Results.sample_window;
            if ~isempty(sample_window)
                mask = ticks >= sample_window(1) & ticks <= sample_window(2);
                ticks = ticks(mask);
            end
            if ~isempty(timeShifts)
                ticks = timeShifts.shiftTimes(ticks);
            end
            if p.Results.timeInSeconds
                ticks = double(ticks) / tsi.fs;
            end
            h = Neuropixel.Utils.rugplot(ticks, 'side', p.Results.side, 'Color', p.Results.Color, 'expand_limits', p.Results.expand_limits);
            set(h, 'XLimInclude', 'off');
        end
        
        function markActiveRegions(tsi, varargin)
            % marks regions that are densely covered with trials (ignoring gaps between trials
            p = inputParser();
            p.addParameter('maxPauseSec', 20, @isscalar); % in samples
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('Color', [0 0 1], @isvector);
            p.addParameter('LineSpecStart', '-', @ischar);
            p.addParameter('LineSpecStop', '--', @ischar);
            p.parse(varargin{:});
            
            [idxStart, idxStop, trialStartStop] = tsi.computeActiveRegions('maxPauseSec', p.Results.maxPauseSec); %#ok<*PROPLC>
            timeShiftSpec = p.Results.time_shifts;
            if ~isempty(timeShiftSpec)
                idxStart = timeShiftSpec.shiftTimes(idxStart);
                idxStop = timeShiftSpec.shiftTimes(idxStop);
            end
            hold on;
            for iT = 1:numel(idxStart)
                str = sprintf('Trial %d-%d', trialStartStop(iT, 1), trialStartStop(iT, 2));
                xline(double(idxStart(iT)) / tsi.fs, p.Results.LineSpecStart, '', 'Color', p.Results.Color, 'LineWidth', 0.5);
                xline(double(idxStop(iT))  / tsi.fs, p.Results.LineSpecStop, str, 'Color', p.Results.Color, 'LineWidth', 0.5);
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
