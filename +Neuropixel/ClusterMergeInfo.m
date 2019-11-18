classdef ClusterMergeInfo < handle
    
    properties(Transient, SetAccess=protected) % don't save KS to disk
        ks
    end
    
    properties(SetAccess=protected) 
        merges (:, 1) cell = {}; % list of sets of cluster_ids, with only one cluster
    end
    
    properties(Dependent)
        new_cluster_ids
        nMerges
    end
    
    methods
        function mi = ClusterMergeInfo(ks)
            assert(isa(ks, 'Neuropixel.KilosortDataset'));
            mi.ks = ks;
        end
        
        function setKilosortDataset(mi, ks)
            assert(isa(ks, 'Neuropixel.KilosortDataset'));
            mi.ks = ks;
        end
        
        function n = get.nMerges(mi)
            n = numel(mi.merges);
        end
        
        function ids = get.new_cluster_ids(mi)
            ids = cellfun(@(clusters) min(clusters), mi.merges);
        end
        
        function postMergeUpdate(mi)
            mi.ensureConsistent();
        end
        
        function clear(mi)
            mi.merges = {};
            mi.postMergeUpdate();
        end
        
        function ensureConsistent(mi)
            maskKeep = true(mi.nMerges, 1);
            for iM = 1:mi.nMerges
                mi.merges{iM} = unique(uint32(mi.merges{iM}(:))');
                maskKeep(iM) = numel(mi.merges{iM}) >= 2;
            end
            mi.merges = mi.merges(maskKeep);
        end
        
        function mergeInds = findMergesForClusters(mi, cluster_ids)
            mask = false(mi.nMerges, 1);
            for m = 1:mi.nMerges
                mask(m) = any(ismember(uint32(cluster_ids), mi.merges{m}));
            end
            mergeInds = find(mask);
        end
        
        function other_cluster_ids = listClustersToBeMergedWith(mi, cluster_ids)
            mergeInds = mi.findMergesForClusters(cluster_ids);
            other_cluster_ids = zeros(0, 1, 'uint32');
            for iM = 1:numel(mergeInds)
                other_cluster_ids = union(other_cluster_ids, mi.merges{mergeInds(iM)});
            end
            other_cluster_ids = setdiff(other_cluster_ids, uint32(cluster_ids));
        end
        
        function mask = isClusterPendingSomeMerge(mi, cluster_ids)
            % returns true for each cluster involved with some merge
            mask = false(numel(cluster_ids), 1);
            for m = 1:mi.nMerges
                mask = mask | ismember(uint32(cluster_ids), mi.merges{m});
            end
        end
        
        function [mask, mergeInds] = isClusterPendingMergeWith(mi, cluster_ids, with_cluster_ids)
            % mask is size(cluster_ids) and true if pending merge with any in with_cluster_ids
            mergeInds = unique(mi.findMergesForClusters(with_cluster_ids));
            mask = false(size(cluster_ids));
            for iM = 1:numel(mergeInds)
                mask = mask | ismember(uint32(cluster_ids), mi.merges{mergeInds(iM)});
            end
        end
        
        function tf = areClustersInSingleMerge(mi, cluster_ids)
            mergeInd = mi.findMergesForClusters(cluster_ids(1));
            if isempty(mergeInd)
                tf = false;
            else
                tf = all(ismember(uint32(cluster_ids), mi.merges{mergeInd}));
            end
        end
        
        function addToMerge(mi, mergeInd, cluster_ids)
            mi.merges{mergeInd} = Neuropixel.Utils.makecol(union(mi.merges{mergeInd}, uint32(cluster_ids)));
            mi.postMergeUpdate();
        end
        
        function mergeMerges(mi, mergeInds, addToNewMerge)
            % merge 2 or more merges by ind
            if islogical(mergeInds)
                mergeInds = find(mergeInds);
            end
            
            newMerge = uint32(addToNewMerge);
            for iM = 1:numel(mergeInds)
                newMerge = union(newMerge, mi.merges{mergeInds(iM)});
            end
            
            keep = setdiff(1:mi.nMerges, mergeInds);
            mi.merges = mi.merges(keep);
            mi.merges{end+1} = Neuropixel.Utils.makecol(unique(newMerge));
            
            mi.postMergeUpdate();
        end
            
        function mergeClusters(mi, cluster_ids)
            mergeInds = mi.findMergesForClusters(cluster_ids);
            if isempty(mergeInds)
                mi.merges{end+1} = Neuropixel.Utils.makecol(unique(uint32(cluster_ids)));
                mi.postMergeUpdate();
            else
                mi.mergeMerges(mergeInds, cluster_ids);
            end
        end
        
        function ensureNotBeingMergedWith(mi, cluster_ids, with_cluster_ids)
            [mask, mergeInds] = mi.isClusterPendingMergeWith(cluster_ids, with_cluster_ids);
            if any(mask)
                for iM = 1:numel(mergeInds)
                    mi.merges{mergeInds(iM)} = setdiff(mi.merges{mergeInds(iM)}, uint32(cluster_ids));
                end
                
                mi.postMergeUpdate();
            end
        end
        
        function excludeFromMerges(mi, cluster_ids)
            maskChange = false(mi.nMerges, 1);
            for m = 1:mi.nMerges
                if any(ismember(cluster_id, mi.merges{m}))
                    mi.merges{m} = setdiff(mi.merges{m}, uint32(cluster_ids));
                    maskChange(m) = true;
                end
            end
            
            if any(maskChange)
                mi.postMergeUpdate();
            end
        end
        
        function removeMerge(mi, mergeInds)
            mi.merges(mergeInds) = [];
            mi.postMergeUpdate();
        end
        
        function [mergeInds, found] = findMergeByNewClusterId(mi, new_cluster_ids)
            [found, mergeInds] = ismember(new_cluster_ids, mi.new_cluster_ids);
        end
        
        function cluster_ids = listClustersWillMergeIntoNewClusterId(mi, new_cluster_ids)
            [mergeInds, found] = mi.findMergeByNewClusterId(new_cluster_ids);
            cluster_ids = mi.listClusterWillMerge(mergeInds(found));
        end
        
        function cluster_ids = listClustersWillMerge(mi, mergeInds)
            cluster_ids = unique(cat(1, mi.merges{mergeInds}));
        end
        
        function cluster_ids_removed_by_merge = listClusterIdsRemovedByApplyingMerges(mi, merge_inds)
            if nargin < 2
                merge_inds = 1:mi.nMerges;
            end
            cluster_ids_removed_by_merge = setdiff(cat(2, mi.merges{merge_inds}), mi.new_cluster_ids);
        end
            
    end
end