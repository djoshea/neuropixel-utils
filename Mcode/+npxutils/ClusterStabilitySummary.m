classdef ClusterStabilitySummary < handle & matlab.mixin.Copyable
    % a data structure for looking at each clusters firing patterns over trials
    % built by KilosortSegmentedDatset .buildClusterStabilitySummary()
    
    properties
        % Trial id of each trial included (nTrials x 1).
        trial_ids(:,1) uint32
        % Whether each trial has data.
        trial_has_data(:,1) logical
        % Condition id of each trial included (numeric or categorical).
        condition_ids(:,1)
        
        % Sample idx of each trial start window, used for x axis in plotting.
        idxStart(:,1) uint64
        % Sampling rate.
        fs(1,1)
        % Time bin width in ms.
        tBinWidth(1,1)
        
        % Id of each cluster (nClusters x 1).
        cluster_ids(:,1) uint32
        
        % nTrials x nClusters x nTime bins.
        rates(:,:)
        
        % Where the file splits occur.
        fileBoundaries(:,1) uint64
        % The names of the split files.
        fileNames(:,1) string
    end
    
    properties(Dependent)
        % Number of trials.
        nTrials
        % Number of clusters.
        nClusters
        % Number of time bins.
        nTimeBins
        % Distinct non-missing condition Ids.
        condition_list
        conditionIndInList
    end
    
    methods
        function n = get.nTrials(this)
            n = size(this.rates, 1);
        end
        
        function n = get.nClusters(this)
            n = size(this.rates, 2);
        end
        
        function n = get.nTimeBins(this)
            n = size(this.rates, 3);
        end
        
        function uc = get.condition_list(this)
            mask = ~ismissing(this.condition_ids);
            uc = unique(this.condition_ids(mask));
        end
        
        function t = getStartTimeSec(this)
            t = double(this.idxStart) / double(this.fs);
        end
        
        function clusterInd = lookup_clusterIds(this, cluster_ids)
            % Get indexes for given cluster IDs.
            %
            % clusterInd = lookup_clusterIds(this, cluster_ids)
            %
            % Cluster_ids (uint32) is a list of cluster IDs to look up.
            %
            % Returns a corresponding array of numeric indexes.
            arguments
                this
                cluster_ids uint32
            end
            [tf, clusterInd] = ismember(cluster_ids, this.cluster_ids);
            if ~all(tf)
                badClusterIds = unique(cluster_ids(~tf));
                error('Clusters not found: %s', mat2str(badClusterIds));
            end
        end
        
        function [conditionInd, condition_list] = lookup_conditionIdsInList(this, condition_ids)
            condition_list = this.condition_list;
            [tf, conditionInd] = ismember(condition_ids, condition_list);
            conditionInd(~tf) = NaN;
        end
        
        function inds = get.conditionIndInList(this)
            inds = this.lookup_conditionIdsInList(this.condition_ids);
        end
        
        function meanFR = computeMeanFRByCondition(this, varargin)
            p = inputParser();
            p.addParameter('ignoreZeroEdges', true, @islogical);
            p.parse(varargin{:});
            
            conditions = this.condition_list;
            meanFR = zeros(numel(conditions), this.nClusters);
            
            for iC = 1:numel(conditions)
                ratesC = this.rates(this.condition_ids == conditions(iC) & this.trial_has_data, :);
                for iClu = 1:this.nClusters
                    if any(ratesC(:, iClu))
                        idxStart = find(ratesC(:, iClu), 1, 'first');
                        idxStop = find(ratesC(:, iClu), 1, 'last');
                        meanFR(iC, iClu) = mean(ratesC(idxStart:idxStop, iClu));
                    end
                end
            end
        end
        
        function plotByCondition(this, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('condition_ids', this.condition_list, @isvector);
            p.addParameter('smoothing', 1, @isscalar);
            p.parse(varargin{:});
            condition_ids = p.Results.condition_ids;
            
            if nargin < 2 || isempty(cluster_ids)
                cluster_ids = this.cluster_ids;
            end
            
            if iscellstr(condition_ids) %#ok<ISCLSTR>
                condition_ids = string(condition_ids);
            end
            clusterInd = this.lookup_clusterIds(cluster_ids);
            
            [conditionIndToPlot, condition_list] = this.lookup_conditionIdsInList(condition_ids); %#ok<*PROPLC>
            conditionIndByTrial = this.conditionIndInList;
            nCondsTotal = numel(condition_list);
            x = this.getStartTimeSec();
            cmap = npxutils.internal.graphics.distinguishable_colors(nCondsTotal);
            
            for iClu = 1:numel(clusterInd)
                for iCond = 1:numel(conditionIndToPlot)
                    trialMask = this.trial_has_data & conditionIndByTrial == conditionIndToPlot(iCond);
                    ratesThis = this.rates(trialMask, clusterInd(iClu));
                    
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