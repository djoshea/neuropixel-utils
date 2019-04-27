classdef ClusterStabilitySummary < handle & matlab.mixin.Copyable
    % a data structure for looking at each clusters firing patterns over trials
    % built by KilosortSegmentedDatset .buildClusterStabilitySummary()
    properties
        % nTrials x 1
        trial_ids(:, 1) uint32 % trial id of each trial included
        trial_has_data(:, 1) logical
        condition_ids(:, 1) % condition id of each trial included (numeric or categorical)
        
        idxStart(:, 1) uint64 % sample idx of each trial start window, used for x axis in plotting
        fs(1,1) % sampling rate
        tBinWidth(1,1) % time bin width in ms
        
        % nClusters x 1
        cluster_ids(:, 1) uint32 % id of each cluster
        
        % nTrials x nClusters x nTime bins
        rates(:, :)
        
        fileBoundaries(:, 1) uint64 % where the file splits occur
        fileNames(:, 1) string % the names of the split files
    end
    
    properties(Dependent)
        nTrials
        nClusters
        nTimeBins
        condition_list
        conditionIndInList
    end
    
    methods
        function n = get.nTrials(s)
            n = size(s.rates, 1);
        end
        
        function n = get.nClusters(s)
            n = size(s.rates, 2);
        end
        
        function n = get.nTimeBins(s)
            n = size(s.rates, 3);
        end
        
        function uc = get.condition_list(s)
            mask = ~ismissing(s.condition_ids);
            uc = unique(s.condition_ids(mask));
        end

        function t = getStartTimeSec(s)
            t = double(s.idxStart) / double(s.fs);
        end
        
        function clusterInd = lookup_clusterIds(s, cluster_ids)
            [tf, clusterInd] = ismember(uint32(cluster_ids), s.cluster_ids);
            assert(all(tf), 'Clusters not found');
        end  
        
        function [conditionInd, condition_list] = lookup_conditionIdsInList(s, condition_ids)
            condition_list = s.condition_list;
            [tf, conditionInd] = ismember(condition_ids, condition_list);
            conditionInd(~tf) = NaN;
        end
        
        function inds = get.conditionIndInList(s)
            inds = s.lookup_conditionIdsInList(s.condition_ids);
        end
        
        function meanFR = computeMeanFRByCondition(s, varargin)
            p = inputParser();
            p.addParameter('ignoreZeroEdges', true, @islogical);
            p.parse(varargin{:});
            
            conditions = s.condition_list;
            meanFR = zeros(numel(conditions), s.nClusters);
            
            for iC = 1:numel(conditions)
                ratesC = s.rates(s.condition_ids == conditions(iC) & s.trial_has_data, :);
                for iClu = 1:s.nClusters
                    if any(ratesC(:, iClu))
                        idxStart = find(ratesC(:, iClu), 1, 'first');
                        idxStop = find(ratesC(:, iClu), 1, 'last');
                        meanFR(iC, iClu) = mean(ratesC(idxStart:idxStop, iClu));
                    end
                end
            end
        end
        
        function plotByCondition(s, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('condition_ids', s.condition_list, @isvector);
            p.addParameter('smoothing', 1, @isscalar);
            p.parse(varargin{:});
            condition_ids = p.Results.condition_ids;
            
            if nargin < 2 || isempty(cluster_ids)
                cluster_ids = s.cluster_ids;
            end
            
            if iscellstr(condition_ids) %#ok<ISCLSTR>
                condition_ids = string(condition_ids);
            end
            clusterInd = s.lookup_clusterIds(cluster_ids);
            
            [conditionIndToPlot, condition_list] = s.lookup_conditionIdsInList(condition_ids); %#ok<*PROPLC>
            conditionIndByTrial = s.conditionIndInList;
            nCondsTotal = numel(condition_list);
            x = s.getStartTimeSec();
            cmap = Neuropixel.Utils.distinguishable_colors(nCondsTotal);
            
            for iClu = 1:numel(clusterInd)
                for iCond = 1:numel(conditionIndToPlot)
                    trialMask = s.trial_has_data & conditionIndByTrial == conditionIndToPlot(iCond);
                    ratesThis = s.rates(trialMask, clusterInd(iClu));
                    
                    if p.Results.smoothing > 1
                        ratesThis = smooth(ratesThis, p.Results.smoothing);
                    end
                    plot(x(trialMask), ratesThis, '.-', 'Color', cmap(conditionIndToPlot(iCond), :));
                    hold on;
                end
            end
            hold off;
        end
    end
end