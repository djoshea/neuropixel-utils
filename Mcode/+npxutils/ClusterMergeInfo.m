classdef ClusterMergeInfo < handle
    
    properties (Transient, SetAccess=protected)
        % The underlying Kilosort data. Transient, so it is not saved to MAT-files.
        ks
    end
    
    properties (SetAccess=protected)
        % List of sets of cluster_ids, with only one cluster.
        merges (:,1) cell = {};
    end
    
    properties (Dependent)
        new_cluster_ids
        nMerges
    end
    
    methods
        function this = ClusterMergeInfo(ks)
            assert(isa(ks, 'npxutils.KilosortDataset'));
            this.ks = ks;
        end
        
        function setKilosortDataset(this, ks)
            assert(isa(ks, 'npxutils.KilosortDataset'));
            this.ks = ks;
        end
        
        function n = get.nMerges(this)
            n = numel(this.merges);
        end
        
        function ids = get.new_cluster_ids(this)
            ids = cellfun(@(clusters) min(clusters), this.merges);
        end
        
        function postMergeUpdate(this)
            this.ensureConsistent();
        end
        
        function clear(this)
            this.merges = {};
            this.postMergeUpdate();
        end
        
        function ensureConsistent(this)
            maskKeep = true(this.nMerges, 1);
            for iM = 1:this.nMerges
                this.merges{iM} = unique(uint32(this.merges{iM}(:))');
                maskKeep(iM) = numel(this.merges{iM}) >= 2;
            end
            this.merges = this.merges(maskKeep);
        end
        
        function mergeInds = findMergesForClusters(this, cluster_ids)
            mask = false(this.nMerges, 1);
            for m = 1:this.nMerges
                mask(m) = any(ismember(uint32(cluster_ids), this.merges{m}));
            end
            mergeInds = find(mask);
        end
        
        function other_cluster_ids = listClustersToBeMergedWith(this, cluster_ids)
            mergeInds = this.findMergesForClusters(cluster_ids);
            other_cluster_ids = zeros(0, 1, 'uint32');
            for iM = 1:numel(mergeInds)
                other_cluster_ids = union(other_cluster_ids, this.merges{mergeInds(iM)});
            end
            other_cluster_ids = setdiff(other_cluster_ids, uint32(cluster_ids));
        end
        
        function mask = isClusterPendingSomeMerge(this, cluster_ids)
            % returns true for each cluster involved with some merge
            mask = false(numel(cluster_ids), 1);
            for m = 1:this.nMerges
                mask = mask | ismember(uint32(cluster_ids), this.merges{m});
            end
        end
        
        function [mask, mergeInds] = isClusterPendingMergeWith(this, cluster_ids, with_cluster_ids)
            % mask is size(cluster_ids) and true if pending merge with any in with_cluster_ids
            mergeInds = unique(this.findMergesForClusters(with_cluster_ids));
            mask = false(size(cluster_ids));
            for iM = 1:numel(mergeInds)
                mask = mask | ismember(uint32(cluster_ids), this.merges{mergeInds(iM)});
            end
        end
        
        function tf = areClustersInSingleMerge(this, cluster_ids)
            mergeInd = this.findMergesForClusters(cluster_ids(1));
            if isempty(mergeInd)
                tf = false;
            else
                tf = all(ismember(uint32(cluster_ids), this.merges{mergeInd}));
            end
        end
        
        function addToMerge(this, mergeInd, cluster_ids)
            this.merges{mergeInd} = npxutils.internal.makecol(union(this.merges{mergeInd}, uint32(cluster_ids)));
            this.postMergeUpdate();
        end
        
        function mergeMerges(this, mergeInds, addToNewMerge)
            % merge 2 or more merges by ind
            if islogical(mergeInds)
                mergeInds = find(mergeInds);
            end
            
            newMerge = uint32(addToNewMerge);
            for iM = 1:numel(mergeInds)
                newMerge = union(newMerge, this.merges{mergeInds(iM)});
            end
            
            keep = setdiff(1:this.nMerges, mergeInds);
            this.merges = this.merges(keep);
            this.merges{end+1} = npxutils.internal.makecol(unique(newMerge));
            
            this.postMergeUpdate();
        end
        
        function mergeClusters(this, cluster_ids)
            mergeInds = this.findMergesForClusters(cluster_ids);
            if isempty(mergeInds)
                this.merges{end+1} = npxutils.internal.makecol(unique(uint32(cluster_ids)));
                this.postMergeUpdate();
            else
                this.mergeMerges(mergeInds, cluster_ids);
            end
        end
        
        function ensureNotBeingMergedWith(this, cluster_ids, with_cluster_ids)
            [mask, mergeInds] = this.isClusterPendingMergeWith(cluster_ids, with_cluster_ids);
            if any(mask)
                for iM = 1:numel(mergeInds)
                    this.merges{mergeInds(iM)} = setdiff(this.merges{mergeInds(iM)}, uint32(cluster_ids));
                end
                
                this.postMergeUpdate();
            end
        end
        
        function excludeFromMerges(this, cluster_ids)
            maskChange = false(this.nMerges, 1);
            for m = 1:this.nMerges
                if any(ismember(cluster_id, this.merges{m}))
                    this.merges{m} = setdiff(this.merges{m}, uint32(cluster_ids));
                    maskChange(m) = true;
                end
            end
            
            if any(maskChange)
                this.postMergeUpdate();
            end
        end
        
        function removeMerge(this, mergeInds)
            this.merges(mergeInds) = [];
            this.postMergeUpdate();
        end
        
        function [mergeInds, found] = findMergeByNewClusterId(this, new_cluster_ids)
            [found, mergeInds] = ismember(new_cluster_ids, this.new_cluster_ids);
        end
        
        function cluster_ids = listClustersWillMergeIntoNewClusterId(this, new_cluster_ids)
            [mergeInds, found] = this.findMergeByNewClusterId(new_cluster_ids);
            cluster_ids = this.listClusterWillMerge(mergeInds(found));
        end
        
        function cluster_ids = listClustersWillMerge(this, mergeInds)
            cluster_ids = unique(cat(1, this.merges{mergeInds}));
        end
        
        function cluster_ids_removed_by_merge = listClusterIdsRemovedByApplyingMerges(this, merge_inds)
            if nargin < 2
                merge_inds = 1:this.nMerges;
            end
            cluster_ids_removed_by_merge = setdiff(cat(2, this.merges{merge_inds}), this.new_cluster_ids);
        end
        
    end
end
