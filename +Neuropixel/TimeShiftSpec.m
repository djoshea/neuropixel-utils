classdef TimeShiftSpec < handle
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
        
        function mat = as_matrix(spec)
            mat = cat(2, spec.idxStart, spec.idxStop, spec.idxShiftStart);
        end
        
        function str = as_string(spec)
            mat = spec.as_matrix();
            str = mat2str(mat);
        end
        
        function dur = get.intervalDurations(spec)
            dur = uint64(int64(spec.idxStop) - int64(spec.idxStart));
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
            srcIdx = spec.computeSourceIndices();
            [~, time_src] = ismember(times, srcIdx);
        end
        
        function constrainToNumSamples(spec, nSamplesSource)
            maskKeep = spec.idxStart <= uint64(nSamplesSource);
            spec.idxStart = spec.idxStart(maskKeep);
            spec.idxStop = min(spec.idxStop(maskKeep), uint64(nSamplesSource));
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
    end
    
    methods(Static)
        function spec = buildToExciseGaps(idxStart, idxStop)
            regionDurations = int64(idxStop) - int64(idxStart);
            regionUpdatedStarts = int64(idxStart);
            regionUpdatedStarts(1) = 1;
            for iR = 2:numel(idxStart)
                regionUpdatedStarts(iR) = regionUpdatedStarts(iR-1) + regionDurations(iR-1);
            end
            
            spec = Neuropixel.TimeShiftSpec(idxStart, idxStop, regionUpdatedStarts);
        end
    end
    
end