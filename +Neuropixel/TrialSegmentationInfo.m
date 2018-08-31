classdef TrialSegmentationInfo
    % simple utility for describing how to segment trials in KiloSort

    properties
        trialId(:, 1) uint32
        conditionId(:, 1) uint32
        idxStart(:, 1) uint64 % in samples
        idxStop(:, 1) uint64 % in samples, may not be used
    end

    methods
        function tsi = TrialSegmentationInfo(nTrials)
            if nargin < 1
                nTrials = 0;
            end

            tsi.trialId = zeros(nTrials, 1, 'uint32');
            tsi.conditionId = zeros(nTrials, 1, 'uint32');
            tsi.idxStart = zeros(nTrials, 1, 'uint64');
            tsi.idxStop = zeros(nTrials, 1, 'uint64');
        end
    end

end
