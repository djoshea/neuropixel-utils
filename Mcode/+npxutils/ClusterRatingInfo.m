classdef ClusterRatingInfo < handle & matlab.mixin.Copyable
    
    properties (Constant, Hidden)
        defaultRatingValueSet = ["unrated", "noise", "lowfr", "unstable", "good"]';
    end
    
    properties (Transient, SetAccess=protected) % don't save KS to disk
        ks
    end
    
    properties (SetAccess=protected)
        cluster_ids (:,1) uint32
        subgroupSet (:,1) string
        
        ratingValueSet (:,1) string;
    end
    
    properties
        % List of ratings corresponding to cluster_ids.
        ratings (:,1) categorical % nClusters
        
        % Is cluster i usable within subgroup j?
        % This is an nClusters x nSubgroups sized array.
        usableWithin (:,:) logical % 
        
        % Is cluster i stable across subgroups j and j+1?
        % This is an nClusters x nSubgroups-1 array indicating whether cluster 
        % is stable across the two adjacent subgroups.
        stableAcross (:,:) logical 
        
        % Include cutoff spikes?
        % An nClusters-long vector.
        includeCutoffSpikes (:, 1) logical
    end
    
    properties (Dependent)
        nClusters
        nSubgroups
        defaultRating
        
        nClustersUnrated
    end
    
    methods
        function this = ClusterRatingInfo(varargin)
            % Construct a new object.
            %
            % obj = npxutils.ClusterRatingInfo(...)
            %
            % Accepts a list of name/value options. Options:
            %
            %   subgroupSet: (string) required
            %   ks:
            %   app:
            %   ratingValueSet: optional
            %
            % One of either ks or app is required.
            p = inputParser();
            % specify this
            p.addParameter('subgroupSet', "All", @isstring);
            
            % specify one of:
            p.addParameter('ks', [], @(x) true);
            p.addParameter('app', [], @(x) isempty(x) || isvector(x));
            
            % optional
            p.addParameter('ratingValueSet', this.defaultRatingValueSet, @isstring);
            
            p.parse(varargin{:})
            
            this.ratingValueSet = p.Results.ratingValueSet;
            
            this.ks = p.Results.ks;
            if isempty(this.ks)
                cluster_ids = p.Results.cluster_ids;
                assert(~isempty(cluster_ids));
            else
                this.ks.load('loadFeatures', false, 'loadBatchwise', false, 'loadCutoff', false);
                cluster_ids = this.ks.cluster_ids;
                assert(~isempty(cluster_ids));
            end
            subgroupSet = p.Results.subgroupSet;
            if isempty(subgroupSet)
                subgroupSet = "";
            end
            this.initializeBlank(cluster_ids, subgroupSet);
        end
        
        function initializeBlank(this, cluster_ids, subgroupSet)
            % Initialize to blank data for a given set of clusters
            %
            % initializeBlank(this, cluster_ids, subgroupSet)
            this.cluster_ids = cluster_ids;
            this.subgroupSet = string(subgroupSet);
            
            this.ratings = repmat(this.defaultRating, this.nClusters, 1);
            this.usableWithin = true(this.nClusters, this.nSubgroups);
            this.stableAcross = true(this.nClusters, this.nSubgroups-1);
            this.includeCutoffSpikes = true(this.nClusters, 1);
        end
        
        function rating = convertToRatingCategorical(this, value)
            % Convert a string rating value to categorical using this' mapping
            %
            % rating = convertToRatingCategorical(this, value)
            %
            % Value (string) is the rating to convert.
            %
            % Returns the ratings as categorical, using this' rating value set.
            % They will be ordinal.
            rating = categorical(string(value), this.ratingValueSet, 'Ordinal', true);
        end
        
        function setKilosortDataset(this, ks)
            % Set the Kilosort dataset
            arguments
                this
                ks npxutils.KilosortDataset
            end
            if ~ks.isLoaded
                ks.load('loadBatchwise', false, 'loadFeatures', false);
            end
            this.ks = ks;
            assert(~isempty(ks.cluster_ids));
            if ~isequal(this.cluster_ids, ks.cluster_ids)
                
                % update my properties to avoid mismatch
                [maskKeep, indSet] = ismember(this.cluster_ids, ks.cluster_ids);
                indSet = indSet(maskKeep);
                
                nDropped = nnz(~maskKeep);
                uniq_ratings_dropped = unique(string(this.ratings(~maskKeep)));
                
                missing = ~ismember(ks.cluster_ids, this.cluster_ids);
                nMissing = nnz(missing);
                warning('KilsortDataset cluster_ids do not match ClusterRatingInfo''s. Missing ratings for %d. Dropping %d with ratings %s', ...
                    nMissing, nDropped, strjoin(uniq_ratings_dropped, '+'));
                
                ratings = this.ratings; %#ok<*PROPLC>
                usableWithin = this.usableWithin;
                stableAcross = this.stableAcross;
                includeCutoffSpikes = this.includeCutoffSpikes;
                this.initializeBlank(ks.cluster_ids, this.subgroupSet);
                
                this.ratings(indSet) = ratings(maskKeep);
                this.usableWithin(indSet, :) = usableWithin(maskKeep, :);
                this.stableAcross(indSet, :) = stableAcross(maskKeep, :);
                this.includeCutoffSpikes(indSet, :) = includeCutoffSpikes(maskKeep, :);
            end
        end
        
        function n = get.nClusters(this)
            n = numel(this.cluster_ids);
        end
        
        function n = get.nSubgroups(this)
            n = max(1, numel(this.subgroupSet));
        end
        
        function rating = get.defaultRating(this)
            rating = this.convertToRatingCategorical(this.defaultRatingValueSet(1));
        end
        
        function [clusterInds, cluster_ids] = lookup_clusterIds(this, cluster_ids)
            % Lookup clusters from a logical mask or numeric list
            %
            % [clusterInds, cluster_ids] = lookup_clusterIds(this, cluster_ids)
            %
            % cluster_ids may be either a logical index into this' cluster_ids
            % list, or a numeric list of cluster ids (not indexes).
            %
            % Returns:
            %   clusterInds - list of corresponding cluster indexes
            %   cluster_ids - numeric cluster ids
            if islogical(cluster_ids)
                cluster_ids = this.cluster_ids(cluster_ids);
            end
            [tf, clusterInds] = ismember(cluster_ids, this.cluster_ids);
            assert(all(tf), 'Some cluster ids were not found in r.cluster_ids');
        end
        
        function [subgroupInds, subgroups] = lookup_subgroups(this, subgroups)
            if isnumeric(subgroups)
                subgroupInds = subgroups;
                subgroups = this.subgroupSet(subgroupInds);
            else
                [tf, subgroupInds] = ismember(string(subgroups), this.subgroupSet);
                assert(all(tf), 'Some subgroups were not found in r.subgroupSet');
                subgroups = this.subgroupSet(subgroupInds);
            end
        end
        
        function [ratings, usableWithin, stableAcross, includeCutoffSpikes] = lookupClusterRatings(this, cluster_ids)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            ratings = this.ratings(cluster_inds);
            
            if nargout > 1
                usableWithin = this.usableWithin(cluster_inds, :);
            end
            if nargout > 2
                stableAcross = this.stableAcross(cluster_inds, :);
            end
            if nargout > 3
                includeCutoffSpikes = this.includeCutoffSpikes(cluster_inds);
            end
        end
        
        function tf = lookupClusterSubgroupUsableWithin(this, cluster_ids, subgroups)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            subgroup_inds = this.lookup_subgroups(subgroups);
            assert(max(subgroup_inds) <= this.nSubgroups);
            
            tf = this.usableWithin(cluster_inds, subgroup_inds);
        end
        
        function tf = lookupClusterSubgroupStableAcross(this, cluster_ids, subgroups)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            subgroup_inds = this.lookup_subgroups(subgroups);
            assert(max(subgroup_inds) < this.nSubgroups);
            
            tf = this.stableAcross(cluster_inds, subgroup_inds);
        end
        
        function setClusterRatings(this, cluster_ids, ratings, usableWithin, stableAcross, includeCutoffSpikes)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            nClu = numel(cluster_inds);
            nSub = this.nSubgroups;
            
            if ~isempty(ratings)
                this.ratings(cluster_inds) = this.convertToRatingCategorical(ratings);
            end
            
            if nargin > 3 && ~isempty(usableWithin)
                if size(usableWithin, 1) == 1
                    usableWithin = repmat(usableWithin, nClu, 1);
                end
                if size(usableWithin, 2) == 1
                    usableWithin = repmat(usableWithin, 1, nSub);
                end
                this.usableWithin(cluster_inds, :) = usableWithin;
            end
            
            if nargin > 4 && ~isempty(stableAcross)
                if size(stableAcross, 1) == 1
                    stableAcross = repmat(stableAcross, nClu, 1);
                end
                if size(stableAcross, 2) == 1
                    stableAcross = repmat(stableAcross, 1, nSub-1);
                end
                this.stableAcross(cluster_inds, :) = stableAcross;
            end
            
            if nargin > 5 && ~isempty(includeCutoffSpikes)
                this.includeCutoffSpikes(cluster_inds) = includeCutoffSpikes;
            end
        end
        
        function tf = lookupClusterIncludeCutoffSpikes(this, cluster_ids)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            tf = this.includeCutoffSpikes(cluster_inds);
        end
        
        function setClusterIncludeCutoffSpikes(this, cluster_ids, includeCutoffSpikes)
            this.setClusterRatings(cluster_ids, [], [], [], includeCutoffSpikes);
        end
        
        function setClusterSubgroupUsableWithin(this, cluster_ids, subgroups, usableWithin)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            subgroup_inds = this.lookup_subgroups(subgroups);
            
            assert(size(usableWithin, 2) == 1 || size(usableWithin, 2) == numel(subgroups));
            this.usableWithin(cluster_inds, subgroup_inds) = usableWithin;
        end
        
        function setClusterSubgroupStableAcross(this, cluster_ids, subgroups, stableAcross)
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            subgroup_inds = this.lookup_subgroups(subgroups);
            assert(all(subgroup_inds < this.nSubgroups));
            
            assert(size(stableAcross, 2) == 1 || size(stableAcross, 2) == numel(subgroups));
            this.stableAcross(cluster_inds, subgroup_inds) = stableAcross;
        end
        
        function countsByRatingSubgroup = computeClusterCounts(this, ratingValueSet, countUsableOnlyMask)
            % Compute cluster counts, possibly filtered by rating value set.
            %
            % countsByRatingSubgroup = computeClusterCounts(this, ratingValueSet, countUsableOnlyMask)
            %
            % countsByRatingSubgroup is numel(ratingValueSet) x nSubgroups count of clusters
            % that have that rating and are listed as usableWithin that subgroup
            
            if nargin < 2 || isempty(ratingValueSet)
                ratingValueSet = this.ratingValueSet;
            end
            ratingValueSet = this.convertToRatingCategorical(ratingValueSet);
            if nargin < 3 || isempty(countUsableOnlyMask)
                ratingsMaskCountAlways = false(this.nClusters, 1);
            else
                ratingsMaskCountAlways = ~ismember(this.ratings, ratingValueSet(countUsableOnlyMask));
            end
            
            nRatingValues = numel(ratingValueSet);
            countsByRatingSubgroup = zeros(nRatingValues, this.nSubgroups);
            
            for iS = 1:this.nSubgroups
                % count a rating if cluster marked usable, or rated something where we count it usable or not
                ratingsCountable = this.usableWithin(:,iS) | ratingsMaskCountAlways;
                countsByRatingSubgroup(:,iS) = histcounts(this.ratings(ratingsCountable), ratingValueSet);
            end
        end
        
        function n = get.nClustersUnrated(this)
            unrated = this.convertToRatingCategorical("unrated");
            n = nnz(this.ratings == unrated);
        end
        
        function [countsByRatingSubgroup, nClustersPostMerge] = ...
                computeClusterCountsAfterApplyingMerges(this, mergeInfo, ratingValueSet, countUsableOnlyMask)
            % Compute cluster counts after applying merges.
            %
            % [countsByRatingSubgroup, nClustersPostMerge] = ...
            %     computeClusterCountsAfterApplyingMerges(this, mergeInfo, ratingValueSet, countUsableOnlyMask)
            %
            % countsByRatingSubgroup is numel(ratingValueSet) x nSubgroups count
            % of clusters that have that rating and are listed as usableWithin
            % that subgroup.
            
            if nargin < 3
                ratingValueSet = [];
            end
            if nargin < 4
                countUsableOnlyMask = [];
            end
            
            % apply merge on copy
            this = copy(this);
            this.apply_cluster_merge(mergeInfo);
            nClustersPostMerge = this.nClusters;
            
            % then count clusters
            countsByRatingSubgroup = this.computeClusterCounts(ratingValueSet, countUsableOnlyMask);
        end
        
        function [nUnrated, nClustersPostMerge, cluster_ids_unrated] = ...
                computeClusterUnratedCountAfterApplyingMerges(this, mergeInfo)
            unrated = this.convertToRatingCategorical("unrated");
            
            this = copy(this);
            this.apply_cluster_merge(mergeInfo);
            nClustersPostMerge = this.nClusters;
            mask_unrated = any(this.ratings == unrated, 2);
            nUnrated = nnz(mask_unrated);
            cluster_ids_unrated = this.cluster_ids(mask_unrated);
        end
        
        function cluster_ids = listClusterIdsUsableWithinSubgroup(this, subgroup, ratingsAccepted)
            % List cluster IDs that are usable within a subgroup.
            %
            % cluster_ids = listClusterIdsUsableWithinSubgroup(r, subgroup, ratingsAccepted)
            %
            % Lists all cluster_ids that are:
            % - rated as one of the ratings in ratingsAccepted
            % - usable within the subgroup listed
            %
            % See also:
            % listClusterIdsUsableAcrossSubgroupsWithRating if aggregating across multiple subgroups
            
            subgroup_ind = sort(this.lookup_subgroups(subgroup), 'ascend');
            mask_rating_accepted = ismember(this.ratings, this.convertToRatingCategorical(ratingsAccepted));
            mask_usable_all = all(this.usableWithin(:, subgroup_ind), 2);
            
            mask = mask_rating_accepted & mask_usable_all;
            cluster_ids = this.cluster_ids(mask);
        end
        
        function [cluster_ids, cluster_ratings] = ...
                listClusterIdsUsableAcrossSubgroupsWithRating(this, subgroups, ratingsAccepted)
            % List cluster IDs usable across subgroups.
            %
            % cluster_ids = listClusterIdsUsableAcrossSubgroupsWithRating(r, subgroups, ratingsAccepted)
            %
            % Lists all cluster_ids that are:
            % - rated as one of the ratings in ratingsAccepted
            % - usable within each subgroup listed,
            % - stable across all interspersed subgroups
            subgroup_inds = sort(this.lookup_subgroups(subgroups), 'ascend');
            mask_rating_accepted = ismember(this.ratings, this.convertToRatingCategorical(ratingsAccepted));
            mask_usable_all = all(this.usableWithin(:, subgroup_inds), 2);
            
            if numel(subgroup_inds) > 1
                mask_stable_across = all(this.stableAcross(:, min(subgroup_inds):max(subgroup_inds)), 2);
            else
                mask_stable_across = true(this.nClusters, 1);
            end
            
            mask = mask_rating_accepted & mask_usable_all & mask_stable_across;
            cluster_ids = this.cluster_ids(mask);
            cluster_ratings = this.ratings(mask);
        end
        
        function apply_cluster_merge(this, mergeInfo)
            ratings = this.ratings; %#ok<*NASGU>
            usableWithin = this.usableWithin;
            stableAcross = this.stableAcross;
            includeCutoffSpikes = this.includeCutoffSpikes;
            
            new_cluster_ids = mergeInfo.new_cluster_ids;
            new_cluster_inds = this.lookup_clusterIds(new_cluster_ids);
            for iM = 1:mergeInfo.nMerges
                cluster_ids = mergeInfo.merges{iM};
                [this_ratings, this_usableWithin, this_stableAcross, this_includeCutoffSpikes] = this.lookupClusterRatings(cluster_ids);
                
                this_new_ind = new_cluster_inds(iM);
                ratings(this_new_ind) = max(this_ratings); % rating is the "max" of merged ratings, which means the best due to the order of defaultRatingValueSet which is used to define the ordering of the categorical
                usableWithin(this_new_ind, :) = any(this_usableWithin, 1); % usable if any is usable
                if ~isempty(stableAcross)
                    stableAcross(this_new_ind, :) = any(this_stableAcross, 1); % stable if any is stable
                end
                includeCutoffSpikes(this_new_ind) = any(this_includeCutoffSpikes); % include cutoff if any includes cutoff
            end
            
            cluster_ids_removed_by_merge = setdiff(cat(2, mergeInfo.merges{:}), new_cluster_ids);
            mask = ~ismember(this.cluster_ids, cluster_ids_removed_by_merge);
            
            % subselect by destination clusters only
            this.ratings = ratings(mask);
            this.usableWithin = usableWithin(mask, :);
            this.stableAcross = stableAcross(mask, :);
            this.includeCutoffSpikes = includeCutoffSpikes(mask);
            this.cluster_ids = this.cluster_ids(mask);
        end
    end
    
    methods (Static)
        function rating = convertToDefaultRatingCategorical(value)
            rating = categorical(string(value), npxutils.ClusterRatingInfo.defaultRatingValueSet, 'Ordinal', true);
        end
        
        function ratings = getVectorUnrated(n)
            ratings = repmat(npxutils.ClusterRatingInfo.convertToDefaultRatingCategorical('unrated'), n, 1);
        end
    end
end
