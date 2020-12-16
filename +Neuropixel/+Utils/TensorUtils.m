classdef TensorUtils
% set of classes for building, manipulating, and computing on
% high-d (or arbitrary-d) matrices easily
%
% Daniel J. O'Shea, 2019

    methods(Static) % Simple internalized utils
        function v = wrapCell(v)
        % v = wrapCell(v)
        % Wrap v as cell {v} if ~iscell(v)
            if ~iscell(v)
                v = {v};
            end
        end

        function varargout = cellvec(N)
            [varargout{1:nargout}] = deal(cell(N, 1));
        end

        function varargout = falsevec(N)
            [varargout{1:nargout}] = deal(false(N, 1));
        end

        function varargout = truevec(N)
            [varargout{1:nargout}] = deal(true(N, 1));
        end

        function varargout = nanvec(N)
            [varargout{1:nargout}] = deal(nan(N, 1));
        end

        function varargout = onesvec(N, varargin)
            [varargout{1:nargout}] = deal(ones(N, 1, varargin{:}));
        end

        function varargout = zerosvec(N, varargin)
            [varargout{1:nargout}] = deal(zeros(N, 1, varargin{:}));
        end

        function vec = makecol( vec )
            % transpose if it's currently a row vector (unless its 0 x 1, keep as is)
            if (size(vec,2) > size(vec, 1) && isvector(vec)) && ~(size(vec, 1) == 0 && size(vec, 2) == 1)
                vec = vec';
            end

            if size(vec, 1) == 1 && size(vec, 2) == 0
                vec = vec';
            end
        end

        function vec = makerow( vec )
        % convert vector to row vector

            % leave size == [1 0] alone too
            if(size(vec, 1) > size(vec,2) && isvector(vec)) && ~(size(vec, 2) == 0 && size(vec, 1) == 1)
                vec = vec';
            end

            if size(vec, 1) == 0 && size(vec, 2) == 1
                vec = vec';
            end
        end

    end


    methods(Static) % Mapping, construction via callback
        function varargout = mapToSizeFromSubs(sz, varargin)
            % t = mapTensor(sz, contentsFn = @(varargin) NaN, asCell = false)
            % build a tensor with size sz by passing subscripts inds to
            % contentsFn(sub1, sub2, ...) maps subscript indices as a vector to the contents
            % asCell == true --> returns cell, asCell == false returns matrix, defaults to false

            p = inputParser;
            p.addRequired('size', @(x) isempty(x) || (isvector(x) && isnumeric(x)));
            p.addOptional('contentsFn', [], @(x) isa(x, 'function_handle'));
            p.addOptional('asCell', false, @islogical);
            p.addOptional('minDims', 2, @isscalar);
            p.parse(sz, varargin{:});
            asCell = p.Results.asCell;
            contentsFn = p.Results.contentsFn;

            if isempty(sz)
                varargout = Neuropixel.Utils.TensorUtils.cellvec(nargout);
                for i = 1:nargout
                    if asCell
                        varargout{i} = {};
                    else
                        varargout{i} = [];
                    end
                end
                return
            end

            sz = Neuropixel.Utils.TensorUtils.expandSizeToNDims(sz, p.Results.minDims);
            nDims = length(sz);
            idxEachDim = arrayfun(@(n) 1:n, sz, 'UniformOutput', false);
            [subsGrids{1:nDims}] = ndgrid(idxEachDim{:});

            if isempty(contentsFn)
                if asCell
                    contentsFn = @(varargin) {};
                else
                    contentsFn = @(varargin) NaN;
                end
            end

            [varargout{1:nargout}] = arrayfun(contentsFn, subsGrids{:}, 'UniformOutput', ~asCell);
        end

        function varargout = map(fn, varargin)
            % works just like cellfun or arrayfun except auto converts each arg
            % to a cell so that cellfun may be used. Returns a cell arrays with
            % the same size as the tensor, although cell arrays containing
            % the same
            for iArg = 1:length(varargin)
                if ~iscell(varargin{iArg})
                    varargin{iArg} = num2cell(varargin{iArg});
                end
            end
            [varargout{1:nargout}] = cellfun(fn, varargin{:}, 'UniformOutput', false);
            % convert scalar numeric cells back to matrices
%             for iArg = 1:numel(varargout);
%                 if all(cellfun(@(x) isnumeric(x) && isscalar(x), varargout{iArg}));
%                     varargout{iArg} = cell2mat(varargout{iArg});
%                 end
%             end
        end

        function results = mapIncludeSubs(fn, varargin)
            % mapWithInds(fn, t1, t2, ...) calls fn(t1(subs), t2(subs), ..., subs) with subs
            % being the subscript indices where the element of t1, t2, etc.
            % was extracted

            for iArg = 1:length(varargin)
                if ~iscell(varargin{iArg})
                    varargin{iArg} = num2cell(varargin{iArg});
                end
            end
            tSubs = Neuropixel.Utils.TensorUtils.containingSubscripts(size(varargin{1}));
            results = cellfun(fn, varargin{:}, tSubs, 'UniformOutput', false);
        end

        function varargout = mapIncludeSubsAndSize(fn, varargin)
            % mapWithInds(fn, t1, t2, ...) calls fn(t1(subs), t2(subs), ..., subs, sz) with subs
            % being the subscript indices where the element of t1, t2, etc.
            % was extracted and sz being size(t1) == size(t2).

            sz = size(varargin{1});
            fnWrap = @(varargin) fn(varargin{:}, sz);
            [varargout{1:nargout}] = Neuropixel.Utils.TensorUtils.mapIncludeSubs(fnWrap, varargin{:});
        end

        function varargout = mapSlices(fn, spanDim, varargin)
            % varargout = mapSlices(fn, spanDims, varargin)
            %
            % this acts like map, calling fn(varargin{1}(ind),varargin{2}(ind))
            % except rather than being called on each element of varargin{:}
            % individually, it is called on slices of the tensor(s) at once. These slices
            % are created by selecting all elements along the dimensions in spanDims and
            % repeating this over each set of subscripts along the other dims.
            % The slices passed to fn will not be squeezed, so they will
            % have singleton dimensions for dim in spanDim. Call
            % .squeezeDims(in, spanDim) to obtain a squeezed slice.
            %
            % The result will be reassembled into a tensor, whose size is determined by
            % the sizes of dimensions not in spanDim. Because the output
            % values will be stored as cell tensor elements, there are no
            % constraints on what these outputs look like

            %sz = size(varargin{1});
            nd = ndims(varargin{1});
            %nArgs = length(varargin);

            % we select individual slices by selecting each along the non-spanned dims
            dim = setdiff(1:nd, spanDim);

            if isempty(dim)
                % when all dims are selected
                tCellArgs = cellfun(@(t) {t}, varargin, 'UniformOutput', false);
            else
                % slice through each of the varargin
                tCellArgs = cellfun(@(t) Neuropixel.Utils.TensorUtils.selectEachAlongDimension(t, dim), ...
                    varargin, 'UniformOutput', false);
            end

            % run the function on each slice
            [resultCell{1:nargout}] = cellfun(fn, tCellArgs{:}, 'UniformOutput', false);

            varargout = resultCell;

            % (old) reassemble the result
            % varargout = cellfun(@(r) Neuropixel.Utils.TensorUtils.reassemble(r, dim), resultCell, 'UniformOutput', false);
        end

        function varargout = mapSlicesInPlace(fn, spanDim, varargin)
            % This function acts very similarly to mapSlices. The only
            % difference is that fn must return outputs that has the same
            % size and shape as it's inputs. Provided this constraint is
            % met, the output will be a tensor that has the same shape as
            % the input tensor. Effectively you loop through the
            % input tensors one slice at a time, transform that slice in
            % place via slice = fn(slice) and rebuild the output tensor one
            % slice at a time.

            [resultCell{1:nargout}] = Neuropixel.Utils.TensorUtils.mapSlices(fn, spanDim, varargin{:});
%             nonSpanDims = Neuropixel.Utils.TensorUtils.otherDims(size(resultCell{1}), spanDim, ndims(varargin{1}));
            
            varargout = cellfun(@cell2mat, resultCell, 'UniformOutput', false);
%             varargout = cellfun(@(r) Neuropixel.Utils.TensorUtils.reassemble(r, nonSpanDims, ndims(varargin{1})), ...
%                 resultCell, 'UniformOutput', false);

%             for iV = 1:numel(varargin)
%                 if ~iscell(varargin{iV})
%                     varargout{iV} = cell2mat(varargout{iV});
%                 end
%             end
        end

        function varargout = mapFromAxisLists(fn, axisLists, varargin)
            % varargout = mapFromAxisLists(fn = @(varargin) ..., axis1, axis2, ...)
            % build a tensor with size sz by setting
            % tensor(i,j,k,...) = fn(axis1{i}, axis2{j}, axis3{k})
            % (or axis(i) for numeric arrays)
            % if fn is omitted returns a cell array of
            % {axis1{i}, axis2{j}, axis3{j}}

            p = inputParser;
            p.addParameter('asCell', true, @islogical);
            p.addParameter('collectArguments', false, @islogical); % if true, call as fn({arg 1, arg 2, arg 3}) rather than fn(arg1, arg2, arg3)
            p.parse(varargin{:});

            if isempty(fn)
                fn = @(varargin) varargin;
            end
            if isa(axisLists, 'function_handle')
                error('Usage: mapFromAxisLists(fn, axisLists, ...)');
            end
            if iscell(axisLists)
                sz = cellfun(@numel, axisLists);
            else
                sz = arrayfun(@numel, axisLists);
            end

            [varargout{1:nargout}] = Neuropixel.Utils.TensorUtils.mapToSizeFromSubs(sz, @indexFn, 'asCell', p.Results.asCell);

            function varargout = indexFn(varargin)
                inputs = Neuropixel.Utils.TensorUtils.cellvec(numel(axisLists));
                for iAx = 1:numel(axisLists)
                    if iscell(axisLists{iAx})
                        inputs{iAx} = axisLists{iAx}{varargin{iAx}};
                    else
                        inputs{iAx} = axisLists{iAx}(varargin{iAx});
                    end
                end

                if p.Results.collectArguments
                    [varargout{1:nargout}] = fn(inputs);
                else
                    [varargout{1:nargout}] = fn(inputs{:});
                end
            end

        end

        function t = buildCombinatorialStructTensor(varargin)
            % Given all struct vector arguments passed in, constructs a
            % tensor struct array containing the merged version of each
            % struct array

            nByArg = cellfun(@numel, varargin);
            fieldsByArg = cellfun(@fieldnames, varargin, 'UniformOutput', false);
            allFields = cat(1, fieldsByArg{:});

            asCells = cell(nargin, 1);
            nArg = numel(varargin);
            for iArg = 1:nArg
                % get as vector with size nByArg(iArg) entries
                flatVector = varargin{iArg}(:);
                % struct2cell returns nFields x nByArg(iArg) cell array,
                % which we transpose
                thisAsCell = struct2cell(flatVector)';
                % which we orient to have size nFields along dim nargin+1
                % and size nByArg(iArg) along dim iArg
                thisAsCellOriented = Neuropixel.Utils.TensorUtils.orientSliceAlongDims(thisAsCell, [iArg, nArg+1]);
                % and expand to be the same size as the full combinatorial tensor
                asCells{iArg} = Neuropixel.Utils.TensorUtils.singletonExpandToSize(thisAsCellOriented, nByArg);
            end

            % combined will be size [nByArg(:) nFieldsTotal)
            combined = cat(nargin+1, asCells{:});
            t = cell2struct(combined, allFields, nArg+1);
        end

        function S = buildCombinatorialStringTensorFromLists(listsByDim, joinWith)
            if nargin < 2
                joinWith = ' ';
            end

            % convert numerics to strings
            for iL = 1:numel(listsByDim)
                if isempty(listsByDim{iL}) || numel(listsByDim{iL}) == 1 && isempty(listsByDim{iL}{1})
                    % just leave as singleton dimension, don't alter string
                    listsByDim{iL} = {''};
                else
                    for iV = 1:numel(listsByDim{iL})
                        if isnumeric(listsByDim{iL}{iV})
                            listsByDim{iL}{iV} = num2str(listsByDim{iL}{iV});
                        end

                        % and add the joinWith string to the end
                        if iL < numel(listsByDim)
                            listsByDim{iL}{iV} = [listsByDim{iL}{iV}, joinWith];
                        end
                    end
                end
            end

            S = Neuropixel.Utils.TensorUtils.mapFromAxisLists(@horzcat, listsByDim);
        end

        function t = mapCatToTensor(fn, varargin)
            % like arrayfun, except with UniformOutput false, and the
            % results will be concatenated together using cellmat to form a
            % tensor
            t = Neuropixel.Utils.TensorUtils.map(fn, varargin{:});
            if iscell(t)
                t = cell2mat(t);
            end
        end
    end

    methods(Static) % Dimensions and sizes
        function out = emptyWithSameType(in, szOut)
            % create out as size szOut with same type as in (either cell or
            % nans)
            if iscell(in)
                out = cell(szOut);
            else
                out = nan(szOut);
            end
        end

        function tf = isvector(vec)
            % like isvector except works for high-d vectors as well
            sz = size(vec);
            tf = nnz(sz ~= 1) == 1 && ~any(sz == 0);
        end

        function sz = sizeMultiDim(t, dims)
            % sz = sizeMultiDim(t, dims) : sz(i) = size(t, dims(i))
            szAll = Neuropixel.Utils.TensorUtils.expandSizeToNDims(size(t), max(dims));
            sz = arrayfun(@(d) szAll(d), dims);
        end

        function sz = sizeOtherDims(t, excludeDims)
            % sz = sizeOtherDims(t, excludeDims)
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(size(t), excludeDims);
            sz = Neuropixel.Utils.TensorUtils.sizeMultiDim(t, otherDims);
        end

        function sz = expandScalarSize(sz)
            % if sz (size) is a scalar, make it into a valid size vector by
            % appending 1 to the end. i.e. 3 --> [3 1]
            if isempty(sz)
                sz = [0 0];
            elseif isscalar(sz)
                sz = [sz 1];
            end
        end

        function sz = sizeNDims(t, nDims)
            sz = Neuropixel.Utils.TensorUtils.expandSizeToNDims(size(t), nDims);
        end

        % pads sz with 1s to make it nDims length
        function sz = expandSizeToNDims(sz, nDims)
            szPad = ones(nDims - numel(sz), 1);
            sz(end+1:nDims) = szPad;
        end

        function tf = compareSizeVectors(sz1, sz2)
            nDims = max(numel(sz1), numel(sz2));
            tf = isequal(Neuropixel.Utils.TensorUtils.expandSizeToNDims(sz1, nDims), ....
                Neuropixel.Utils.TensorUtils.expandSizeToNDims(sz2, nDims));
        end

        function other = otherDims(sz, dims, ndims)
            % otherDims(t, dims, [ndims]) returns a list of dims in t NOT in dims
            % e.g. if ndims(t) == 3, dims = 2, other = [1 3]
            if nargin < 3
                ndims = max(numel(sz), max(dims));
            end
            allDims = 1:ndims;
            other = Neuropixel.Utils.TensorUtils.makerow(setdiff(allDims, dims));
        end

        function d = firstNonSingletonDim(t)
            sz = size(t);
            d = find(sz > 1, 1, 'first');
            if isempty(d)
                d = 1;
            end
        end
    end

    methods(Static) % Mask generation and mask utilities
       function idx = vectorMaskToIndices(mask)
            if islogical(mask)
                idx = Neuropixel.Utils.TensorUtils.makecol(find(mask));
            else
                idx = Neuropixel.Utils.TensorUtils.makecol(mask);
            end
        end

        function mask = vectorIndicesToMask(idx, len)
            if islogical(idx)
                mask = Neuropixel.Utils.TensorUtils.makecol(idx);
            else
                if ~isscalar(len) && isvector(len), len = numel(len); end % assume its a mask of the same size
                mask = Neuropixel.Utils.TensorUtils.falsevec(len);
                mask(idx) = true;
            end
        end

        function maskNew = subselectVectorMask(maskOrig, maskAnd)
            % given logical maskOrig
            existingKeep = find(maskOrig);
            assert(numel(maskAnd) == numel(existingKeep), 'Subselection mask does not match nnz of original mask');
            keep = Neuropixel.Utils.TensorUtils.vectorIndicesToMask(existingKeep(maskAnd), numel(maskOrig));
            maskNew = maskOrig;
            maskNew(~keep) = false;
        end

        function inflated = inflateMaskedTensor(maskedTensor, dims, masks, fillWith)
            % inflated = inflateMaskedTensor(maskedTensor, dims, masks, fillWith=NaN)
            % takes maskedTensor, which has been formed by slicing some
            % original tensor by selecting with mask along each dimension
            % in dims, and returns the original tensor where the masked out
            % rows, cols, etc. are filled with fill (defaults to NaN)
            % if dims is a scalar, masks is a logical or numeric vector to use
            % for selecting along dim. If dim is a vector, select is a cell array of
            % vectors to be used for selecting along dim(i)

            assert(~islogical(dims), 'Arguments out of order');
            if ~iscell(masks)
                masks = {masks};
            end
            for i = 1:numel(masks)
                if ~islogical(masks{i})
                    masks{i} = Neuropixel.Utils.TensorUtils.vectorIndicesToMask(masks{i}, max(masks{i}));
                end
            end

            assert(all(cellfun(@(x) islogical(x) && isvector(x), masks)), ...
                'Mask must be (cell of) logical vectors. Use vectorIndicesToMask to convert.');

            if numel(masks) == 1
                masks = repmat(masks, numel(dims), 1);
            end
            assert(numel(masks) == numel(dims), 'Number of dimensions must match number of masks provided');
            masks = Neuropixel.Utils.TensorUtils.makecol(masks);

            nInflatedVec = cellfun(@numel, masks);
            nSelectedVec = cellfun(@nnz, masks);

            % check size along selected dims
            assert(isequal(Neuropixel.Utils.TensorUtils.makecol(Neuropixel.Utils.TensorUtils.sizeMultiDim(maskedTensor, dims)), nSelectedVec), ...
                'Size along each masked dimension must match nnz(mask)');

            % compute inflated size
            inflatedSize = size(maskedTensor);
            inflatedSize(dims) = nInflatedVec;

            if nargin < 4
                if iscell(maskedTensor)
                    fillWith = {[]};
                elseif isstruct(maskedTensor)
                    flds = fieldnames(maskedTensor);
                    args = cell(numel(flds)*2, 1);
                    args(1:2:end) = flds;
                    fillWith = struct(args{:});
                elseif isstring(maskedTensor)
                    fillWith = string(missing);
                else
                    fillWith = NaN;
                end
            end
            if ischar(fillWith)
                fillWith = {fillWith};
            end
            inflated = repmat(fillWith, inflatedSize);

            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectAlongDimension(inflatedSize, dims, masks);
            inflated(maskByDim{:}) = maskedTensor;
        end

        function idxFull = indicesIntoMaskToOriginalIndices(idxIntoMasked, mask)
            maskInds = find(mask);
            idxFull = maskInds(idxIntoMasked);
        end

    end

    methods(Static) % Indices and subscripts

        function t = containingLinearInds(sz)
            % build a tensor with size sz where each element contains the linear
            % index it would be accessed at, e.g. t(i) = i
            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);
            t = Neuropixel.Utils.TensorUtils.mapToSizeFromSubs(sz, @(varargin) sub2ind(sz, varargin{:}), false);
        end

        function t = containingSubscripts(sz, asCell)
            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);

            % asCell == true means each element is itself a cell rather then a vector of
            % subscripts
            if nargin < 2
                asCell = false;
            end

            % build a tensor with size sz where each element contains the vector
            % of subscripts it would be accessed at, e.g. t(i) = i
            if asCell
                t = Neuropixel.Utils.TensorUtils.mapToSizeFromSubs(sz, @(varargin) varargin, true);
            else
                t = Neuropixel.Utils.TensorUtils.mapToSizeFromSubs(sz, @(varargin) [varargin{:}]', true);
            end
        end

        function subsCell = ndgridCell(sz)
            % works like ndgrid except captures all outputs

            args = arrayfun(@(s) 1:s, sz, 'UniformOutput', false);
            [subsCell{1:ndims(sz)}] = ndgrid(args{:});
        end

        function t = containingSubscriptsCatAlongFirstDimension(sz)
            % t is a tensor whose size is [numel(sz) sz]
            % t(dim, i, j, k) is the subscript of i,j,k on dimension dim.
            % e.g. if dim == 2, t(2, i, j, k, ...) == j
            t = cell2mat(shiftdim(Neuropixel.Utils.TensorUtils.containingSubscripts(sz), -1));
        end

        function mat = ind2subAsMat(sz, inds)
            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);

            % sz is the size of the tensor
            % mat is length(inds) x length(sz) where each row contains ind2sub(sz, inds(i))

            ndims = length(sz);
            subsCell = cell(ndims, 1);

            [subsCell{:}] = ind2sub(sz, Neuropixel.Utils.TensorUtils.makecol(inds));

            mat = [subsCell{:}];
        end

        function inds = subMat2Ind(sz, mat)
            % inds = subMat2Ind(sz, mat)
            % sz is the size of the tensor
            % mat is length(inds) x length(sz) where each row contains a set of subscripts
            % as would be returned by ind2sub(sz, inds(i))
            % subMat2Ind essentially converts back to linear indices using
            % sub2ind and is the inverse of ind2subAsMat

            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);

            ndims = length(sz);
            % DO NOT UNCOMMENT. THIS WILL BREAK THINGS SINCE SZ CAN HAVE
            % SZ(1) == 1.
             if ndims == 2 && sz(2) == 1 && size(mat, 2) == 1
                 ndims = 1;
             end
            subsCell = arrayfun(@(dim) mat(:, dim), 1:ndims, 'UniformOutput', false);

            inds = sub2ind(sz, subsCell{:});
        end
    end

    methods(Static) % Selection Mask generation
        function maskByDim = maskByDimCell(sz)
            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);

            % get a cell array of selectors into each dim that would select
            % every element if used via t(maskByDim{:})
            maskByDim = arrayfun(@(n) true(n, 1), sz, 'UniformOutput', false);
        end

        % the next few methods accept a dim and select argument
        % if dim is a scalar, select is a logical or numeric vector to use
        % for selecting along dim. If dim is a vector, select is a cell array of
        % vectors to be used for selecting along dim(i)
        function maskByDim = maskByDimCellSelectAlongDimension(sz, dim, select)
            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);

            % get a cell array of selectors into each dim that effectively select
            % select{i} along dim(i). These could be used by indexing a tensor t
            % via t(maskByDim{:}) --> se selectAlongDimension
            if ~iscell(select)
                select = {select};
            end

            assert(length(dim) == length(select), 'Number of dimensions must match length of select mask cell array');
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCell(sz);
            maskByDim(dim) = select;
        end

        function mask = maskSelectAlongDimension(sz, dim, select)
            sz = Neuropixel.Utils.TensorUtils.expandScalarSize(sz);

            % return a logical mask where for tensor with size sz
            % we select t(:, :, select, :, :) where select acts along dimension dim

            mask = false(sz);
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectAlongDimension(sz, dim, select);
            mask(maskByDim{:}) = true;
        end

        function maskByDim = maskByDimCellSelectPartialFromOrigin(sz, dim, select)
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCell(sz);
            maskByDim(dim) = arrayfun(@(n) 1:n, select, 'UniformOutput', false);
        end

        function mask = maskSelectPartialFromOrigin(sz, dim, select)
            % select the first select(i) items from dim(i)
            % for i = 1:numel(sz)
            selectCell = arrayfun(@(n) 1:n, select, 'UniformOutput', false);
            mask = Neuropixel.Utils.TensorUtils.maskSelectAlongDimension(sz, dim, selectCell);
        end
    end

    methods(Static) % Squeezing along particular dimensions
        function tsq = squeezeDims(t, dims)
            % like squeeze, except only collapses singleton dimensions in list dims
            siz = size(t);
            dims = dims(dims <= ndims(t));
            dims = dims(siz(dims) == 1);
            siz(dims) = []; % Remove singleton dimensions.
            siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
            tsq = reshape(t,siz);
        end

        function newDimIdx = shiftDimsPostSqueeze(sz, squeezeDims, dimsToShift)
            assert(isvector(sz), 'First arg must be size');
            % when squeezing along squeezeDims, the positions of dims in
            % dimsToShift will change. The new dim idx will be returned
            origDims = 1:length(sz);
            squeezeDims = squeezeDims(squeezeDims <= length(sz));
            remainDims = setdiff(origDims, squeezeDims);

            [~, newDimIdx] = ismember(dimsToShift, remainDims);
            newDimIdx(newDimIdx==0) = NaN;
        end

        function tsq = squeezeOtherDims(t, dims)
            other = Neuropixel.Utils.TensorUtils.otherDims(size(t), dims);
            tsq = Neuropixel.Utils.TensorUtils.squeezeDims(t, other);
        end
    end

    methods(Static) % Regrouping, Nesting, Selecting, Reshaping
        function tCell = regroupAlongDimension(t, dims)
            % tCell = regroupAlongDimension(t, dims)
            % returns a cell tensor of tensors, where the outer tensor is over
            % the dims in dims. Each inner tensor is formed by selecting over
            % the dims not in dims.
            %
            % e.g. if size(t) = [nA nB nC nD] and dims is [1 2],
            % size(tCell) = [nA nB] and size(tCell{iA, iB}) = [nC nD]

            tCell = Neuropixel.Utils.TensorUtils.squeezeSelectEachAlongDimension(t, dims);
            tCell = Neuropixel.Utils.TensorUtils.squeezeOtherDims(tCell, dims);
        end

        function tCell = nestedRegroupAlongDimension(t, dimSets)
            % this forms cells of cells of cells of ... of tensors
            assert(iscell(dimSets), 'dimSets must be a cell array of dimension sets');

            dimSets = Neuropixel.Utils.TensorUtils.makecol(cellfun(@Neuropixel.Utils.TensorUtils.makecol, dimSets, 'UniformOutput', false));
            allDims = cell2mat(dimSets);

            assert(length(unique(allDims)) == length(allDims), ...
                'A dimension was included in multiple dimension sets');

            otherDims = Neuropixel.Utils.TensorUtils.otherDims(size(t), allDims);
            if ~isempty(otherDims)
                dimSets{end+1} = otherDims;
            end

            tCell = inner(t, dimSets);
            return;

            function tCell = inner(t, dimSets)
                if length(dimSets) == 1
                    % special case, no grouping, just permute dimensions and
                    % force to be cell
                    tCell = permute(t, dimSets{1});
                    if ~iscell(tCell)
                        tCell = num2cell(tCell);
                    end
                elseif length(dimSets) == 2
                    % last step in recursion, call final regroup
                    tCell = Neuropixel.Utils.TensorUtils.regroupAlongDimension(t, dimSets{1});
                else
                    % call inner on each slice of dimSets{1}
                    remainingDims = Neuropixel.Utils.TensorUtils.otherDims(size(t), dimSets{1});
                    tCell = Neuropixel.Utils.TensorUtils.mapSlices(@(t) mapFn(t, dimSets), ...
                        remainingDims, t);
                    tCell = Neuropixel.Utils.TensorUtils.squeezeOtherDims(tCell, dimSets{1});
                end
            end

            function tCell = mapFn(t, dimSets)
                remainingDimSets = cellfun(...
                    @(dims) Neuropixel.Utils.TensorUtils.shiftDimsPostSqueeze(size(t), dimSets{1}, dims), ...
                    dimSets(2:end), 'UniformOutput', false);
                tCell = inner(Neuropixel.Utils.TensorUtils.squeezeDims(t, dimSets{1}), remainingDimSets);
            end

        end

        function [res, maskByDim] = selectAlongDimension(t, dim, select, squeezeResult)
            if nargin < 4
                squeezeResult = false;
            end
            sz = size(t);
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectAlongDimension(sz, dim, select);
            res = t(maskByDim{:});

            if squeezeResult
                % selectively squeeze along dim
                res = Neuropixel.Utils.TensorUtils.squeezeDims(res, dim);
            end
        end

        function t = assignIntoTensorAlongDimension(t, assignThis, dims, select)
            sz = size(t);
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectAlongDimension(sz, dims, select);
            t(maskByDim{:}) = assignThis;
        end

        function selected = selectMaskAlongMultipleDimensions(t, dims, mask)
            % mask logical mask or vector of linear inds into dimensions in
            % dims. selected will be size along other dims x nnz(mask)
            sz = size(t);
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(sz, dims);
            t_reshape = Neuropixel.Utils.TensorUtils.reshapeByConcatenatingDims(t, {otherDims, dims});
            selected = t_reshape(:, mask);

            selected = reshape(selected, [sz(otherDims) size(selected, 2)]);
        end

        function t = assignIntoTensorAlongMultipleDimensionsByMask(t, assignThis, dims, mask)
            % mask logical mask or vector of linear inds into dimensions in
            % dims. assignThis should scalar or be size(t)
            sz = size(t);
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(sz, dims);
            t_reshape = Neuropixel.Utils.TensorUtils.reshapeByConcatenatingDims(t, {otherDims, dims});

            if isscalar(assignThis)
                t_reshape(:, mask) = assignThis;
            else
                assignThis = Neuropixel.Utils.TensorUtils.reshapeByConcatenatingDims(assignThis, {otherDims, dims});
                t_reshape(:, mask) = assignThis;
            end

            t = Neuropixel.Utils.TensorUtils.undoReshapeByConcatenatingDims(t_reshape, {otherDims, dims}, sz);
        end

        function sel = selectAlongDimensionWithNaNs(t, dim, select, varargin)
            assert(numel(dim) == 1, 'Must be single dimension');
            nanMask = isnan(select(:));
            selNoNan = Neuropixel.Utils.TensorUtils.selectAlongDimension(t, dim, select(~nanMask), varargin{:});
            sel = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(selNoNan, dim, ~nanMask);
        end

        function [res, maskByDim] = squeezeSelectAlongDimension(t, dim, select)
            % select ind along dimension dim and squeeze() the result
            % e.g. squeeze(t(:, :, ... ind, ...))

            [res, maskByDim] = Neuropixel.Utils.TensorUtils.selectAlongDimension(t, dim, select, true);
        end

        function tCell = selectEachAlongDimension(t, dim, squeezeEach)
            % returns a cell array tCell such that tCell{i} = selectAlongDimension(t, dim, i)
            % optionally calls squeeze on each element
            if nargin < 3
                squeezeEach = false;
            end

            sz = size(t);

            % generate masks by dimension that are equivalent to ':'
            %maskByDimCell = Neuropixel.Utils.TensorUtils.maskByDimCell(sz);

            dimMask = true(ndims(t), 1);
            dimMask(dim) = false;
            szResult = sz;
            szResult(dimMask) = 1;

            % oh so clever
            tCell = Neuropixel.Utils.TensorUtils.mapToSizeFromSubs(szResult, 'minDims', max(dim), 'asCell', true, ...
                'contentsFn', @(varargin) Neuropixel.Utils.TensorUtils.selectAlongDimension(t, dim, ...
                varargin(dim), squeezeEach));
        end

        function tCell = squeezeSelectEachAlongDimension(t, dim)
            % returns a cell array tCell such that tCell{i} = squeezeSelectAlongDimension(t, dim, i)
            tCell = Neuropixel.Utils.TensorUtils.selectEachAlongDimension(t, dim, true);
        end

        function out = splitAlongDimensionBySubscripts(t, dim, outSz, subs)
            % given an n-dimensional tensor t, split t into a cell with size outsz
            % each cell at location (s1, s2, ...) will contain a piece of t selected along dimension dim taking all 
            % slices t(..., i, ...) where subs(i, :) == (s1, s2, ...). Sort of like accumarray except concatenating over slices
            % where any subs(i, :) is NaN or 0, this slice of t will be discarded
            
            assert(size(subs, 1) == size(t, dim), 'Number of rows of subs must match size of t along dim');
            assert(size(subs, 2) == numel(outSz), 'Number of elements in outSz should match number of columns in subs');
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(size(t), dim);
            szT = size(t);
            
            mask = ~any(isnan(subs) | subs == 0, 2);
            if any(~mask)
                t = Neuropixel.Utils.TensorUtils.selectAlongDimension(t, dim, mask);
                subs = subs(mask, :);
            end
            
            C = size(subs, 2);
            subsExp = nan(numel(t), C);
            for iC = 1:C
                colOrient = Neuropixel.Utils.TensorUtils.orientSliceAlongDims(subs(:, iC), dim);
                subcolExp = Neuropixel.Utils.TensorUtils.repmatAlongDims(colOrient, otherDims, szT(otherDims));
                subsExp(:, iC) = subcolExp(:);
            end
            
            outSz = Neuropixel.Utils.TensorUtils.expandScalarSize(outSz);
            
            szInnerArgs = num2cell(szT);
            szInnerArgs{dim} = [];
            
            % must be sorted to preserve order in accumarray
            if size(subsExp, 2) > 1
                [subsSorted,subSortIdx] = sortrows(subsExp,[2,1]);
            else
                [subsSorted,subSortIdx] = sort(subsExp);
            end
            if isempty(subsSorted)
                % accumarray won't return a cell if it never evaluates the function handle
                out = repmat({zeros(szInnerArgs{:}, 'like', t)}, outSz);
            else
                out = accumarray(subsSorted, t(subSortIdx), outSz, @(x) {x}, {zeros(0, 1, 'like', t)});
                % reshape out{:} back to tensors 
                out = cellfun(@(x) reshape(x, szInnerArgs{:}), out, 'UniformOutput', false);
            end
        end
        
        function out = splitAlongDimensionByIndex(t, dim, which, outSz)
            % see splitBySubscripts which handles subscripts for which instead of linear indices
            % and is dramatically faster using accumarray
            
            if nargin < 4
                outSz = max(which);
            end
            
            assert(numel(which) == size(t, dim), 'Number of rows of subs must match size of t along dim');
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(size(t), dim);
            szT = size(t);
            
            mask = ~(isnan(which) | which == 0);
            if any(~mask)
                t = Neuropixel.Utils.TensorUtils.selectAlongDimension(t, dim, mask);
                subs = which(mask, :);
            end
            
            colOrient = Neuropixel.Utils.TensorUtils.orientSliceAlongDims(which, dim);
            subcolExp = Neuropixel.Utils.TensorUtils.repmatAlongDims(colOrient, otherDims, szT(otherDims));
            subsExp = subcolExp(:);

            outSz = Neuropixel.Utils.TensorUtils.expandScalarSize(outSz);
            
            % must be sorted to preserve order in accumarray
            [subsSorted,subSortIdx] = sort(subsExp);
            out = accumarray(subsSorted, t(subSortIdx), outSz, @(x) {x}, zeros(0, 1, 'like', t));
            
            % reshape out{:} back to tensors 
            szInnerArgs = num2cell(szT);
            szInnerArgs{dim} = [];
            out = cellfun(@(x) reshape(x, szInnerArgs{:}), out, 'UniformOutput', false);
        end
        
        % it is unclear to me now what this function did, that is not already done by cell2mat.
%         function t = reassemble(tCell, dim, nd)
%             % given a tCell in the form returned by selectEachAlongDimension
%             % return the original tensor. So if tCell is size 3 x 1 x 1 x 6 and tCell{:} is 1 x 4 x 5 x 1,
%             % the output will be 3 x 4 x 5 x 6 with element i,j,k,l drawn from tCell{i,1,1,l}(1,j,k,1)
%             % dim is the dimensions that tCell spans, the other dimensions are spanned by each individual tCell{i}
%             
% 
%             %             if nargin < 3
%             %                 nd = ndims(tCell);
%             %             end
%             szOuter = size(tCell);
%             szOuter = [szOuter ones(1, nd - length(szOuter))];
%             szInner = size(tCell{1});
%             szInner = [szInner ones(1, nd - length(szInner))];
% 
%             % dimMask(i) true i
%             dimMask = false(nd, 1);
%             dimMask(dim) = true;
% 
%             % compute size of result t
%             % use outerDims when its in dim, innerDims when it isn't
%             szT = nan(1, ndims(tCell));
%             szT(dimMask) = szOuter(dimMask);
%             szT(~dimMask) = szInner(~dimMask);
% 
%             % rebuild t by grabbing the appropriate element from tCell
%             %subs = Neuropixel.Utils.TensorUtils.containingSubscripts(szT);
%             t = Neuropixel.Utils.TensorUtils.mapToSizeFromSubs(szT, @getElementT, true);
% 
%             function el = getElementT(varargin)
%                 [innerSubs, outerSubs] = deal(varargin);
%                 % index with dim into tt, non-dim into tt{i}
%                 [outerSubs{~dimMask}] = deal(1);
%                 [innerSubs{dimMask}] = deal(1);
%                 tEl = tCell{outerSubs{:}};
%                 if iscell(tEl)
%                     el = tEl{innerSubs{:}};
%                 else
%                     el = tEl(innerSubs{:});
%                 end
%             end
%         end

        function vec = flatten(t)
            vec = Neuropixel.Utils.TensorUtils.makecol(t(:));
        end

        function mat = flattenAlongDimension(t, dim)
            % returns a 2d matrix where mat(i, :) is the flattened vector of tensor
            % values from each t(..., i, ...) where i is along dim

            nAlong = size(t, dim);
            nWithin = numel(t) / nAlong;
            if iscell(t)
                mat = cell(nAlong, nWithin);
            else
                mat = nan(nAlong, nWithin);
            end

            sqMask = Neuropixel.Utils.TensorUtils.maskByDimCell(size(t));
            for iAlong = 1:nAlong
                sqMask{dim} = iAlong;
                within = t(sqMask{:});
                mat(iAlong, :) = within(:);
            end
        end

        function tCell = flattenAlongDimensionAsCell(t, dim)
            % returns a cell array of length size(t, dim)
            % where each element is the flattened vector of tensor
            % values from each t(..., i, ...) where i is along dim
            tCell = Neuropixel.Utils.TensorUtils.regroupAlongDimension(t, dim);
            for iAlong = 1:length(tCell)
                tCell{iAlong} = tCell{iAlong}(:);
            end
        end

        function [out, which] = catWhich(dim, varargin)
            % works like cat, but returns a vector indicating which of the
            % inputs each element of out came from
            out = cat(dim, varargin{:});
            if nargout > 1
                which = cell2mat(Neuropixel.Utils.TensorUtils.makecol(cellfun(@(in, idx) idx*ones(size(in, dim), 1), varargin, ...
                    num2cell(1:numel(varargin)), 'UniformOutput', false)));
            end
        end

        function paddedCell = expandToSameSizeAlongDims(dims, varargin)
            nd = max(max(dims), max(cellfun(@ndims, varargin)));
            szMat = cell2mat(cellfun(@(x) Neuropixel.Utils.TensorUtils.sizeNDims(x, nd), Neuropixel.Utils.TensorUtils.makecol(varargin), 'UniformOutput', false));

            szPadTo = max(szMat, [], 1);
            padFn = @(x) Neuropixel.Utils.TensorUtils.expandOrTruncateToSize(x, dims, szPadTo(dims));

            paddedCell = cellfun(padFn, Neuropixel.Utils.TensorUtils.makecol(varargin), 'UniformOutput', false);
        end
        
        function paddedCell = expandToSameSize(varargin)
            nd = max(cellfun(@ndims, varargin));
            dims = 1:nd;
            paddedCell = Neuropixel.Utils.TensorUtils.expandToSameSizeAlongDims(dims, varargin{:});
        end
        
        function out = catPad(dim, varargin)
            % works like cat, except pads each element to be the same size
            % along dimension dim before concatenating

            nd = max(cellfun(@ndims, varargin));
            otherDims = setdiff(1:nd, dim);
            paddedCell = Neuropixel.Utils.TensorUtils.expandToSameSizeAlongDims(otherDims, varargin{:});
            out = cat(dim, paddedCell{:});
        end

        function [out, tvec] = concatenateAlignedToCommonTimeVector(timeCell, dataCell, timeDim, catDim, varargin)
            p = inputParser();
            p.addParameter('timeVec', [], @isvector);
            p.addParameter('timeDelta', [], @isscalar);
            p.parse(varargin{:});

            if isempty(p.Results.timeVec)
                timeDelta = p.Results.timeDelta;
                if isempty(timeDelta)
                    timeDelta = nanmedian(timeCell{1});
                end

                % first find common time vector
                tMin = min(cellfun(@min, timeCell));
                tMax = min(cellfun(@max, timeCell));

                tvec = tMin:timeDelta:tMax;
            else
                % use provided time vector
                tvec = p.Results.timeVec;
            end

            % expand / align each dataCell to live on time vector
            nd = max(max([timeDim, catDim]), max(cellfun(@ndims, dataCell)));
            szMat = cell2mat(cellfun(@(x) Neuropixel.Utils.TensorUtils.sizeNDims(x, nd), dataCell, 'UniformOutput', false));

            finalSize = szMat(1, :);
            finalSize(catDim) = sum(szMat(:, catDim));
            finalSize(timeDim) = numel(tvec);
            out = nan(finalSize, 'like', dataCell{1});

            % assume uniform sampling
            Ncat = numel(dataCell);
            for i = 1:Ncat
                ti = timeCell{i};
                idxStart = Neuropixel.Utils.TensorUtils.argMin(abs(ti(1) - tvec));
                idxEnd = idxStart + numel(ti) - 1;

                dims = [catDim, timeDim];
                masks = {i, idxStart:idxEnd};

                out = Neuropixel.Utils.TensorUtils.assignIntoTensorAlongDimension(out, dataCell{i}, dims, masks);
            end
        end

        function dataNew = sliceOrExpandToAlignTimeVector(tvecCurrent, dataCurrent, tvecNew, timeDim)
            assert(size(dataCurrent, timeDim) == numel(tvecCurrent));

            % first slice off ends if needed
            mask = tvecCurrent >= min(tvecNew) & tvecCurrent <= max(tvecNew);
            dataSlice = Neuropixel.Utils.TensorUtils.selectAlongDimension(dataCurrent, timeDim, mask);
            tvecCurrent = tvecCurrent(mask);

            % then expand as neede
            sz = Neuropixel.Utils.TensorUtils.sizeNDims(dataSlice, timeDim);
            sz(timeDim) = numel(tvecNew);
            dataNew = nan(sz, 'like', dataSlice);

            idxStart = Neuropixel.Utils.TensorUtils.argMin(abs(min(tvecCurrent) - tvecNew));
            idxEnd = idxStart + numel(tvecCurrent) - 1;
            dataNew = Neuropixel.Utils.TensorUtils.assignIntoTensorAlongDimension(dataNew, dataSlice, timeDim, idxStart:idxEnd);
        end

        function [out, which] = catWhichIgnoreEmpty(dim, varargin)
            % works like cat, but returns a vector indicating which of the
            % inputs each element of out came from

            isEmpty = cellfun(@isempty, varargin);
            out = cat(dim, varargin{~isEmpty});
            if isempty(out)
                which = [];
                return;
            end
            if nargout > 1
                whichMasked = cell2mat(Neuropixel.Utils.TensorUtils.makecol(cellfun(@(in, idx) idx*ones(size(in, dim), 1), varargin(~isEmpty), ...
                    num2cell(1:nnz(~isEmpty)), 'UniformOutput', false)));

                % whichMasked indexes into masked varargin, reset these to
                % index into the original varargin
                which = Neuropixel.Utils.TensorUtils.indicesIntoMaskToOriginalIndices(whichMasked, ~isEmpty);
            end
        end

        function [out, which] = catInnerDimOverOuterDim(t, innerDim, outerDim)
            % for each entry in t along dimension outerDim, concatenates
            % these entries along dimension, innerDim

            if nargout == 2
                [out, whichCell]  = Neuropixel.Utils.TensorUtils.mapSlices(@(slice) Neuropixel.Utils.TensorUtils.catWhich(innerDim, slice{:}), outerDim, t);
                which = whichCell{1};
            else
                out = Neuropixel.Utils.TensorUtils.mapSlices(@(slice) cat(innerDim, slice{:}), outerDim, t);
            end
        end

        function [out, labelsByDimOut] = reshapeByConcatenatingDims(in, whichDims, labelsByDim)
            % reshapes a tensor by concatenating dims.
            % whichDims is a cell vector indicating which dimensions of in to
            % capture along each dimension of out.
            %
            % For example, suppose in has size s1 x s2 x s3 x s4 x s5.
            % And whichDims = { [ 1 2 ], [3 5], 4 }
            % Then out will be 3-dimensional with size:
            %    (s1*s2) x (s3*s5) x (s4)
            % such that in(i1, i2, i3, i4, i5) will end up at
            %   out(o1, o2, o3), where:
            %       o1 = sub2ind([s1 s2], i1, i2);
            %       o2 = sub2ind([s3 s5], i3, i5);
            %       o3 = sub2ind(s4, i4) == s4;
            % This ensures that along out dim 1, all elements of in-dim 1
            % will remain together for dim2=1, then all elements along
            % in-dim 1 for dim2=2, etc. That is, the first dimension will
            % be swept with the later dimensions constant. This is the
            % reverse of how you would implement this using nested for loops,
            % i.e. [1 2] means for j = 1:size(in, 2), for i = 1:size(in, 1)
            % and is ultimately because Matlab uses column major matrices
            %
            % LabelsByDimOut allows the user to keep track of what pieces
            % of in end up at which positions in out. LabelsByDimIn
            % represent a set of labels along each of the "axes" or
            % dimensions of in. LabelsByDimOut represent a set of labels
            % along each of the 'axes" of dimension out. Those labels are
            % taken to match the labels of in and are concatenated in the
            % same fashion.
            %
            % labelsByDim (optional) is a ndims(in) cell vector, each
            % element i contains a vector or matrix with size(in, i) rows and
            % an arbitrary number of colums. If labelsByDim{i} is a vector,
            % it MUST be a column vector.
            %
            % labelsByDimOut will have a similar format: a ndims(out) cell
            % where each element i contains a matrix with size(out, i)
            % rows. The contents are formed by concatenating the columns
            % taken from each of the labelsByDim{...} cells.
            %
            % The format of labelsByDimIn and labelsByDimOut is that the
            % output can be used as a input to this function again.

            if ~iscell(whichDims)
                whichDims = {whichDims};
            end
            whichDims = Neuropixel.Utils.TensorUtils.makecol(whichDims);

            allDims = cellfun(@(x) x(:), whichDims, 'UniformOutput', false);
            allDims = cat(1, allDims{:});

            szIn = Neuropixel.Utils.TensorUtils.sizeNDims(in, max(allDims));
            ndimsIn = max(max(allDims), ndims(in));

            assert(all(ismember(1:numel(allDims), allDims)), ...
                'whichDims must contain each dim in 1:length(whichDims)');
            % add any trailing dimensions which are missing from the list
            % automatically
            if max(allDims) < ndimsIn
                whichDims = cat(1, whichDims, num2cell(max(allDims)+1:ndimsIn));
                % recompute allDims in case trailing dims were added
                allDims = cellfun(@(x) x(:), whichDims, 'UniformOutput', false);
                allDims = cat(1, allDims{:});
            end

            ndimsOut = length(whichDims);

            szOut = nan(1, ndimsOut);
            for iDimOut = 1:ndimsOut
                szOut(iDimOut) = prod(szIn(whichDims{iDimOut}));
            end

            % order dimensions of out according to whichDims
            out = reshape(permute(in, allDims), szOut);

            % build labels for output dimensions
            if nargout > 1
                szInExpand = Neuropixel.Utils.TensorUtils.expandScalarSize(szIn);
                labelsByDimTemplate = arrayfun(@(dim) 1:szInExpand(dim), 1:numel(szInExpand), 'UniformOutput', false);
                if ~exist('labelsByDim', 'var')
                    labelsByDim = labelsByDimTemplate;
                else
                    % extend labelsByDim to the full length ndims if necessary
                    labelsByDimTemplate(1:numel(labelsByDim)) = labelsByDim;
                    labelsByDim = labelsByDimTemplate;
                end

                szOut = size(out);
                labelsByDimOut = Neuropixel.Utils.TensorUtils.cellvec(ndimsOut);
                for iDimOut = 1:ndimsOut
                    dimsFromIn = whichDims{iDimOut};
                    % temporary storage of the columns of
                    % labelsByDimOut{iDimOut}, before concatenation
                    labelsThisOut = Neuropixel.Utils.TensorUtils.cellvec(numel(dimsFromIn));
                    subsCell = Neuropixel.Utils.TensorUtils.cellvec(numel(dimsFromIn));

                    % get subscripts for each element in the dimsFromIn slice
                    [subsCell{1:numel(dimsFromIn)}] = ind2sub(szIn(dimsFromIn), 1:szOut(iDimOut));
                    for iDimFromIn = 1:numel(dimsFromIn)
                        % check whether it's a correctly sized row vector and convert to
                        % column vector
                        labelsThisIn = labelsByDim{dimsFromIn(iDimFromIn)};
                        if isvector(labelsThisIn) && length(labelsThisIn) == size(in, dimsFromIn(iDimFromIn))
                            labelsThisIn = Neuropixel.Utils.TensorUtils.makecol(labelsThisIn);
                        end
                        % and pull the correct labels by selecting the
                        % correct rows (typically multiple times)
                        labelsThisOut{iDimFromIn} = labelsThisIn(subsCell{iDimFromIn}, :);
                    end
                    % then concatenate these columns together
                    labelsByDimOut{iDimOut} = cat(2, labelsThisOut{:});
                end
            end

        end

        function out = undoReshapeByConcatenatingDims(in, whichDims, szOrig)
            % if is passed in {[1 2]} and orig has 3 dims, whichDims would
            % become {[1 2] 3}

            ndimsOrig = numel(szOrig);

            % process whichDims as in reshapeByConcatenatingDims
            if ~iscell(whichDims)
                whichDims = {whichDims};
            end
            whichDims = Neuropixel.Utils.TensorUtils.makecol(whichDims);
            allDims = cellfun(@(x) x(:), whichDims, 'UniformOutput', false);
            allDims = cat(1, allDims{:});
            assert(all(ismember(1:numel(allDims), allDims)), ...
                'whichDims must contain each dim in 1:length(whichDims)');
            % add any trailing dimensions which are missing from the list
            % automatically
            if max(allDims) < ndimsOrig
                whichDims = cat(1, whichDims, num2cell(max(allDims)+1:ndimsOrig));
                % recompute allDims in case trailing dims were added
                allDims = cellfun(@(x) x(:), whichDims, 'UniformOutput', false);
                allDims = cat(1, allDims{:});
            end

            % intermediate size to "unconcatenate" the dims that were
            % concatenated
            if numel(szOrig) < max(allDims)
                szOrig = [szOrig ones(1, max(allDims)-numel(szOrig))];
            end
            szExpand = szOrig(allDims(1:numel(szOrig)));
            expand = reshape(in, szExpand);

            % reorder dimensions
            out = ipermute(expand, allDims);
        end

        function [out, newDims] = reshapeDimsInPlace(in, whichDims, newSizeInThoseDims)
            % reshapes a tensor by taking the consecutive dimension in
            % whichDims and making these a new size, e.g. tensorizes a
            % flattened dim or set of dims
            sz = size(in);
            assert(isequal(Neuropixel.Utils.TensorUtils.makecol(whichDims), Neuropixel.Utils.TensorUtils.makecol(min(whichDims):max(whichDims))), 'Dims must be consecutive');

            nElements = prod(Neuropixel.Utils.TensorUtils.sizeMultiDim(in, whichDims));
            if nnz(isnan(newSizeInThoseDims)) == 1
                mask = isnan(newSizeInThoseDims);
                fillIn = nElements / prod(newSizeInThoseDims(~mask));
                assert(fillIn == round(fillIn), 'Size specification cannot be satisfied');
                newSizeInThoseDims(mask) = fillIn;
            else
                assert(prod(newSizeInThoseDims) == nElements, 'Size must not change during reshape');
            end
            mind = min(whichDims);
            maxd = max(whichDims);

            szNew = [sz(1:mind-1), newSizeInThoseDims, sz(maxd+1:end)];
            out = reshape(in, szNew);
            newDims = (mind : mind+numel(newSizeInThoseDims)-1)';
        end

        function out = flattenDimsInPlace(in, whichDims)
            newLen = prod(Neuropixel.Utils.TensorUtils.sizeMultiDim(in, whichDims));
            out = Neuropixel.Utils.TensorUtils.reshapeDimsInPlace(in, whichDims, newLen);
        end

        function out = selectSetsAlongDimension(in, dim, selectIdxCell)
            % given a tensor with size s1 x s2 x s3, where lets say dim is 2,
            % selects multiple elements along dimension 2 for every
            % position in dimensions 1, 3. if selectIdxCell has length
            % S2, returns cell with size s1 x S2 x s3

            otherDims = Neuropixel.Utils.TensorUtils.otherDims(size(in), dim);
            each = cellfun(@(idx) Neuropixel.Utils.TensorUtils.squeezeSelectEachAlongDimension(...
                Neuropixel.Utils.TensorUtils.selectAlongDimension(in, dim, idx), otherDims), selectIdxCell, ...
                'UniformOutput', false);
            out = cat(dim, each{:});
        end

        function out = selectSpecificIndicesAlongDimensionEachPosition(in, dim, idxThatDim)
            % typically we do a selection along dim e.g. where dim = 2, using in(:, idx, :)
            % this function enables idx to vary with the position along the
            % other dimensions. If size(in) = in szIn, and
            % idxForEachOtherDim has szIdx, which matches szIn except on
            % dimension dim.
            % for example, if in = [1 2; 3 4], dim = 1, and
            % idxForEachOtherDim = [1 2], out is [1 4]

            maskValid = ~isnan(idxThatDim);
            idxThatDim(~maskValid) = 1;

            function out = selectIdx(slice, idx)
                if isnan(idx)
                    out = NaN;
                else
                    out = slice(idx);
                end
            end
            out = Neuropixel.Utils.TensorUtils.mapSlices(@selectIdx, dim, in, idxThatDim);
            if ~iscell(in)
                out = cell2mat(out);
                out(~maskValid) = NaN;
            else
                if dim == 1
                    out = cat(2, out{:});
                else
                    out = cat(1, out{:});
                end
            end
        end

    end

    methods(Static) % Slice orienting and repmat
        function vec = orientVectorAlongDim(vec, dim)
            vec = shiftdim(Neuropixel.Utils.TensorUtils.makecol(vec), -dim+1);
        end

        % A slice is a selected region of a tensor, in which each dimension
        % either selects all of the elements (via :), or 1 of the elements.
        % This is the 2-d matrix equivalent of a row (slice with spanDim = 1)
        % or a column (slice with spanDim = 2). For a 3-d tensor, a slice
        % along spanDim = [2 3] would look like t(1, :, :);
        function out = orientSliceAlongDims(slice, spanDim)
            % given a slice with D dimensions, orient slice so that it's
            % dimension(i) becomes dimension spanDim(i)

            %szSlice = size(slice);
            ndimsSlice = length(spanDim);
            if ndimsSlice == 1
                % when spanDim is scalar, expect column vector or Neuropixel.Utils.TensorUtils.makecol
                slice = Neuropixel.Utils.TensorUtils.makecol(squeeze(slice));
            end

            ndimsOut = max(2, max(spanDim));
            sliceHigherDims = ndimsSlice+1:ndimsOut;
            nonSpanDims = setdiff(1:ndimsOut, spanDim);

            permuteOrder = nan(1, ndimsOut);
            permuteOrder(spanDim) = 1:ndimsSlice;
            permuteOrder(nonSpanDims) = sliceHigherDims;

            out = permute(slice, permuteOrder);
        end

        function out = orientSliceAlongSameDimsAs(slice, refSlice)
            % out = orientSliceAlongSameDimsAs(slice, refSlice)
            % orient slice such that it has the same shape as refSlice
            spanDim = find(size(refSlice) > 1);
            if ~isempty(spanDim);
                out = Neuropixel.Utils.TensorUtils.orientSliceAlongDims(slice, spanDim);
            else
                out = slice;
            end
        end

        function out = repmatSliceAlongDims(slice, szOut, spanDim)
            % given a slice with ndims(slice) == numel(spanDim),
            % orient that slice along the dimensions in spanDim and repmat
            % it so that the output is size szOut.

            %nSpan = 1:length(spanDim);
            %ndimsOut = length(szOut);

            sliceOrient = Neuropixel.Utils.TensorUtils.orientSliceAlongDims(slice, spanDim);

            repCounts = szOut;
            repCounts(spanDim) = 1;

            out = repmat(sliceOrient, repCounts);
        end

        function t = unshiftdim(t, nshift, ndimsOrig)
            % t = unshiftdim(t, nshift, ndimsOrig)
             if nshift > 0
                 % shiftdim did permute(t,[n+1:ndims(x),1:n])
                 t = ipermute(t, [nshift+1:ndimsOrig,1:nshift]);
             elseif nshift < 0
                 % shiftdim just did a reshape, so we can just shiftdim
                 % back
                 t = shiftdim(t, -nshift);
             end
        end

        function [t, ndimOrig] = shiftdimToFirstDim(t, dim)
            ndimOrig = ndims(t);
            if dim ~= 1
                t = shiftdim(t, dim-1);
            end
        end

        function t = unshiftdimToFirstDim(t, dim, ndimsOrig)
            t = Neuropixel.Utils.TensorUtils.unshiftdim(t, dim-1, ndimsOrig);
        end
    end

    methods(Static) % Size expansion / truncation
        function out = scalarExpandToSize(in, szOut)
            if isscalar(in)
                out = repmat(in, szOut);
            else
                % check that size is okay
                szIn = size(in);
                szOut = Neuropixel.Utils.TensorUtils.expandScalarSize(szOut);
                assert(isequal(szOut, szIn));
                out = in;
            end
        end

        function out = singletonExpandToSize(in, szOut)
            % expand singleton dimensions to match szOut using repmat
            szIn = size(in);
            repCounts = szOut;
            repCounts(szIn ~= 1) = 1;
            out = repmat(in, repCounts);
        end

        function out = expandAlongDims(in, dims, by, fillWith)
            % expand by adding nans or empty cells for cell array
            sz = size(in);
            sz(dims) = sz(dims) + by;

            if nargin < 4 || isempty(fillWith)
                out = Neuropixel.Utils.TensorUtils.emptyWithSameType(in, sz);
            else
                out = repmat(fillWith, sz);
            end

            % build the mask over out to assign in into
            szInDims = Neuropixel.Utils.TensorUtils.sizeMultiDim(in, dims);
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectPartialFromOrigin(size(out), dims, szInDims);

            out(maskByDim{:}) = in;
        end

        function out = repmatAlongDims(in, dims, toSz)
            sz = size(in);
            repmatArg = sz;
            repmatArg(:) = 1;
            repmatArg(dims) = toSz;
            out = repmat(in, repmatArg);
        end

        function t = expandOrTruncateToSize(t, dims, makeSize, fillWith)
            % along each dims(i), truncate or expand with NaN to be size sz(i)

            sz = Neuropixel.Utils.TensorUtils.sizeMultiDim(t, dims);
            makeSize = Neuropixel.Utils.TensorUtils.makerow(makeSize);
            expandBy = makeSize - sz;

            if nargin < 4
                fillWith = [];
            end

            tooSmall = expandBy > 0;
            if any(tooSmall)
                t = Neuropixel.Utils.TensorUtils.expandAlongDims(t, dims(tooSmall), expandBy(tooSmall), fillWith);
            end

            tooLarge = expandBy < 0;
            if any(tooLarge)
                maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectPartialFromOrigin(size(t), dims(tooLarge), makeSize(tooLarge));
                t = t(maskByDim{:});
            end
        end
    end

    methods(Static) % List resampling and shuffling along dimension
        function seedRandStream(seed)
            if nargin == 0
                seed = 'Shuffle';
            end
            s = RandStream('mt19937ar', 'Seed', seed);
            RandStream.setGlobalStream(s);
        end

        % Each of these methods apply to cell tensors containing vector
        % lists of values, or matrices of values where each row is kept together.
        % They will move values (or rows) around among the
        % different cells, i.e. from one list to another, but generally
        % they each preserve the total (row) count within each list.
        % They also typically act along one or more dimensions, in that
        % they will move values among the cells that lie along a particular
        % dimension(s), but not among cells that lie at different positions
        % on the other dimensions. For example, if t is a 2-d matrix,
        % shuffling along dimension 1 would move values among cells that
        % lie along the same row, but not across rows. Here, the row would
        % be referred to as a slice of the matrix.
        %
        % Inputs are not required to have lists be column vectors, but the
        % output will have column vectors in each cell.
        %
        % The various functions differ in how they move the values around
        % and whether or not they sample with replacement, which would
        % allow some values to be replicated multiple times and others to
        % not appear in the final lists.

        function t = listHorzCatContents(varargin)
            % given a set of lists (cells of vectors or matrices), horzcat each lists  contents along the rows
            t = cell(size(varargin{1}));
            for iC = 1:numel(varargin{1})
                elEach = cellfun(@(x) x{iC}, varargin, 'UniformOutput', false);
                t{iC} = cat(2, elEach{:});
            end
        end

        function varargout = listSplitContents(t)
            % given a set of lists (cells of vectors or matrices), horzcat each lists  contents along the rows

            if isempty(t)
                % provide the expected number of outputs, each empty
                nO = nargout;
                varargout = cell(nO, 1);
                [varargout{:}] = deal(cell(size(t)));
                return;
            end

            nO = size(t{1}, 2);
            varargout = cell(nO, 1);
            [varargout{:}] = deal(cell(size(t)));

            for iO = 1:nO
                for iT = 1:numel(t)
                    varargout{iO}{iT} = t{iT}(:, iO);
                end
            end
        end

        function [t] = listShuffleAlongDimension(t, iA, replace)
            % t is a cell tensor of vectors or matrices.
            % For each slice of t along dimension iA, i.e.
            % t(i, j, ..., :, ..., k) where the : is in position iA,
            % shuffle the elements among all cells of that slice,
            % preserving the total number in each cell. If replace is true
            % shuffles with replacement.

            if nargin < 3
                replace = false;
            end

            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@shuffleFn, iA, t);

            function sNew = shuffleFn(s)
                list = Neuropixel.Utils.TensorUtils.combineListFromCells(s);
                if size(list, 1) > 1
                    % randsample with scalar first input thinks we mean
                    % 1:scalar as the population rather than just [scalar]
                    list = list(randsample(size(list, 1), size(list, 1), replace), :, :, :, :);
                end
                sNew = Neuropixel.Utils.TensorUtils.splitListIntoCells(list, cellfun(@numel, s));
            end
        end

        function t = listResampleFromSame(t)
            % resample from replacement within each list, i.e. don't move
            % anything between lists, just take each list and resample with
            % replacement from itself.
            t = Neuropixel.Utils.TensorUtils.map(@(list) list(randsample(size(list, 1), size(list, 1), true), :, :, :, :), t);
        end

        function t = listResampleFromSpecifiedAlongDimension(t, from, iA, replace)
            % resample from replacement from other lists, according to a
            % specified set of cell elements along each slice. from{i} is a
            % specifies the linear index (or indices) of which cells to
            % sample from when building the list for the cell at position i
            % along that slice. If replace is true, samples with
            % replacement.
            %
            % Example: listResampleFromSpecified(t, {1, [1 2]}, 1) where t
            % is a 3 x 2 matrix. For each row i of t, the new list at cell
            % t{i, 1} will be built by sampling from list in t{i, 1}.
            % The new list at cell t{i, 2} will be built by sampling from
            % both t{i, 1} and t{i, 2}.

            if nargin < 4
                replace = true;
            end

            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@resampleFn, iA, t);

            function sNew = resampleFn(s)
                sNew = cell(size(s));
                for i = 1:numel(s)
                    list = Neuropixel.Utils.TensorUtils.combineListFromCells(s(from{i}));
                    sNew{i} = list(randsample(size(list, 1), size(s{i}, 1), replace), :, :, :, :);
                end
            end
        end

        function list = combineListFromCells(t)
            % gathers all lists from t{:} into one long vector/matrix/tensor along the first dimension
            % make empty vectors column vectors to allow cell2mat to work
            [t{cellfun(@isempty, t)}] = deal(nan(0, 1));
            list = cell2mat(t(:));
        end

        function t = splitListIntoCells(list, nPerCell)
            % list is a vector list, nPerCell is a numeric tensor
            % specifying how many elements to place in each list.
            % This undoes combineListFromCells
            if isempty(list)
                % all cells will be empty, special case
                t = Neuropixel.Utils.TensorUtils.cellvec(numel(nPerCell));
            else
                szel = size(list);
                t = mat2cell(list, nPerCell(:), szel(2:end));
            end
            t = reshape(t, size(nPerCell));

        end

        function c = splitAlongDimension(t, dim, nPerCell, squeezeDims)
           % c = splitAlongDimension(t, dim, nPerCell)
           % wrapper around mat2cell that splits t along dimension t so
           % that size(c) == [numel(nPerCell) 1] and size(c{i}, dim) == nPerCell(i)
           % if dim is a vector, then nPerCell may be a cell vector which
           % will independently split along each dimension in dim. If
           % nPerCell is omitted, will be ones(size(t, dim))

           args = Neuropixel.Utils.TensorUtils.cellvec(ndims(t));
           if nargin < 3 || isempty(nPerCell)
               nPerCell = arrayfun(@(dim) Neuropixel.Utils.TensorUtils.onesvec(size(t, dim)), dim, 'UniformOutput', false);
           end
           if nargin < 4
               squeezeDims = false;
           end
           if ~iscell(nPerCell)
               nPerCell = {nPerCell};
           end

           assert(numel(nPerCell) == numel(dim), 'nPerCell must be vector if dim is scalar or cell with one element per dim');
           for iDim = 1:ndims(t)
               [tf, which] = ismember(iDim, dim);
               if tf
                   if isscalar(nPerCell{which}) && nPerCell{which} ~= size(t, iDim)
                       nPerCell_this = repmat(nPerCell{which}, floor(size(t, iDim) / nPerCell{which}), 1);
                       rem = mod(size(t, iDim), nPerCell{which});
                       if rem > 0
                           nPerCell_this(end+1) = rem; %#ok<AGROW>
                       end
                       args{iDim} = nPerCell_this;
                   else
                       args{iDim} = nPerCell{which};
                   end
               else
                   args{iDim} = size(t, iDim);
               end
           end

           c = mat2cell(t, args{:});

           if squeezeDims
               c = cellfun(@(x) Neuropixel.Utils.TensorUtils.squeezeDims(x, dim), c, 'UniformOutput', false);
           end
        end

        function c = splitIntoBlocks(t, dims, nPerDim, squeezeDims)
            if nargin < 4
                squeezeDims = false;
            end

            args = Neuropixel.Utils.TensorUtils.cellvec(ndims(t));

            assert(numel(nPerDim) == numel(dims), 'nPerCell must be vector if dim is scalar or cell with one element per dim');
            for iDim = 1:ndims(t)
                [tf, which] = ismember(iDim, dims);
                s = size(t, iDim);
                if tf
                    nthis = nPerDim(which);
                    nreps = floor(s / nthis);
                    last = rem(s, nthis);
                    if last == 0
                        args{iDim} = repmat(nthis, nreps, 1);
                    else
                        args{iDim} = [repmat(nthis, nreps, 1); last];
                    end
                else
                    args{iDim} = s;
                end
            end

            c = mat2cell(t, args{:});

            if squeezeDims
                c = cellfun(@(x) Neuropixel.Utils.TensorUtils.squeezeDims(x, dim), c, 'UniformOutput', false);
            end
        end


        function t = blockfun(t, blockSz, fn, asTensor)
            % simpler version of matlab's blockproc that works in n-d
            if nargin < 4
                asTensor = false;
            end

            t = Neuropixel.Utils.TensorUtils.splitIntoBlocks(t, 1:numel(blockSz), blockSz);
            t = cellfun(fn, t, 'UniformOutput', false);
            if asTensor
                t = cell2mat(t);
            end
        end

        function t = blockMean(t, blockSz)
            t = Neuropixel.Utils.TensorUtils.blockfun(t, blockSz, @(x) mean(x(:), 'omitnan'), true);
        end

        

    end

    methods(Static) % Multi-dim extensions of any, all, find etc.
        function t = applyFunctionComposedOverSuccessiveDimensions(t, fn, dims)
            for iD = 1:numel(dims)
                t = fn(t, dims(iD));
            end
        end

        function t = anyMultiDim(t, dims)
            % works like any, except operates on multiple dimensions
            for iD = 1:numel(dims)
                t = any(t, dims(iD));
            end
        end

        function t = allMultiDim(t, dims)
            % works like any, except operates on multiple dimensions
            for iD = 1:numel(dims)
                t = all(t, dims(iD));
            end
        end

        function t = nanmaxMultiDim(t, dims)
            % works like any, except operates on multiple dimensions
            for iD = 1:numel(dims)
                t = nanmax(t, [], dims(iD));
            end
        end

        function t = nanminMultiDim(t, dims)
            % works like any, except operates on multiple dimensions
            for iD = 1:numel(dims)
                t = nanmin(t, [], dims(iD));
            end
        end

        function r = nanrangeMultiDim(t, dims)
            r = Neuropixel.Utils.TensorUtils.nanmaxMultiDim(t, dims) - Neuropixel.Utils.TensorUtils.nanminMultiDim(t, dims);
        end

        function t = quantileMultiDim(t, p, dim, placeAlongDim)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            if nargin < 4
                placeAlongDim = dim(1);
            end
            % quantile returns row vectors, so we orient them along placeAlongDim
            shift = -placeAlongDim + 2;
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) shiftdim(quantile(slice(:), p), shift), dim, t);
        end

        function idxTensor = findNAlongDim(t, dim, N, direction)
            % idxTensor = findNAlongDim(t, dim, N, direction)
            % finds the first/last N non-zero values in t along dimension t
            % proceeds in direction 'first' (default) or 'last'
            % if less than N non-zero values are found, idxTensor will be
            % NaN. idxTensor is the same size as t along all dimensions
            % except dim, where it has size N.

            if nargin < 4
                direction = 'first';
            end
            if ~islogical(t)
                t = t ~= 0;
            end

            % create a template for findInner to use with the correct size
            % and orientation
            sizeSlice = Neuropixel.Utils.TensorUtils.onesvec(max(dim, ndims(t)))';
            sizeSlice(dim) = N;
            emptySlice = nan(sizeSlice);

            idxTensor = cell2mat(Neuropixel.Utils.TensorUtils.mapSlices(@findInner, dim, t));

            function w = findInner(v)
                w = emptySlice;
                idx = find(v, N, direction);
                w(1:numel(idx)) = idx;
            end
        end
    end

    methods(Static) % Statistics
        function t = meanMultiDim(t, dims)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) mean(slice(:)), dims, t);
        end

        function t = nanmeanMultiDim(t, dims)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) nanmean(slice(:)), dims, t);
        end

        function t = nansumMultiDim(t, dims)
            for iD = 1:numel(dims)
                t = nansum(t, dims(iD));
            end
        end

        function t = varMultiDim(t, normOpt, dims)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) var(slice(:), normOpt), dims, t);
        end

        function t = nanvarMultiDim(t, normOpt, dims)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) var(slice(:), normOpt, 'omitnan'), dims, t);
        end

        function t = stdMultiDim(t, dims, varargin)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) std(slice(:), varargin{:}, 'omitnan'), dims, t);
        end

        function t = nanstdMultiDim(t, dims, varargin)
            % e.g. if t has size [s1, s2, s3, s4], then  mean(t, [2 3])
            % will compute the mean in slices along dims 2 and 3. the
            % result will have size s1 x 1 x 1 x s4
            t = Neuropixel.Utils.TensorUtils.mapSlicesInPlace(@(slice) std(slice(:), varargin{:}, 'omitnan'), dims, t);
        end

        function t = zscoreMultiDim(t, alongDims)
            % for each subscript in dimension(s) alongDims, computes the mean
            % along all other dimensions and subtracts it. this ensures
            % that the mean along any slice in alongDims will have zero
            % mean. Then normalizes by the std along all other dimensions.
            t = Neuropixel.Utils.TensorUtils.centerSlicesSpanningDimension(t, alongDims);
            stdTensor = Neuropixel.Utils.TensorUtils.nanstdMultiDim(t, alongDims);
            t = bsxfun(@rdivide, t, stdTensor);
        end

        function t = zscoreSoftMultiDim(t, alongDims, denomOffset)
            % for each subscript in dimension(s) alongDims, computes the mean
            % along all other dimensions and subtracts it. this ensures
            % that the mean along any slice in alongDims will have zero
            % mean. Then normalizes by the std along all other dimensions.
            t = Neuropixel.Utils.TensorUtils.centerSlicesSpanningDimension(t, alongDims);
            stdTensor = Neuropixel.Utils.TensorUtils.nanstdMultiDim(t, alongDims) + denomOffset;
            t = bsxfun(@rdivide, t, stdTensor);
        end

        function [y, group] = buildAnovanInputs(valueCell, namesAlongAxes)
            % converts from tensor array of vectors into format needed for
            % Matlab's anovan
            %
            % valueCell is a cell array with multiple axes, one for each factor
            %
            % namesAlongAxes is nAxes x 1 each containing the names of levels along
            % each axis of value cell

            nAx = ndims(valueCell);
            if nAx == 2 && size(valueCell, 2) == 1
                nAx = 1;
            end
            sz = size(valueCell);

            if nargin < 2 || isempty(namesAlongAxes)
                for a = 1:nAx
                    namesAlongAxes{a} = (1:sz(a))';
                end
            end
            assert(numel(namesAlongAxes) == nAx, 'Number of axes must match numel(namesAlongAxes)');

            group = Neuropixel.Utils.TensorUtils.cellvec(nAx);
            for a = 1:nAx
                group{a} = cell(sz);
            end

            % construct the cell of factor level names
            subs = Neuropixel.Utils.TensorUtils.cellvec(nAx);
            for i = 1:numel(valueCell)
                [subs{1:nAx}] = ind2sub(sz, i);
                n = numel(valueCell{i});
                for a = 1:nAx
                    group{a}{i} = repmat(namesAlongAxes{a}(subs{a}), n, 1);
                end

                valueCell{i} = Neuropixel.Utils.TensorUtils.makecol(valueCell{i});
            end

            % collapse over treatments
            y = cat(1, valueCell{:});
            for a = 1:nAx
                group{a} = cat(1, group{a}{:});
            end
        end

        function [idxAlongDim, val] = argMin(t, dim)
            % returns the position along dim of the minimum of each
            % position
            if nargin < 2
                dim = Neuropixel.Utils.TensorUtils.firstNonSingletonDim(t);
            end
            [val, idxAlongDim] = Neuropixel.Utils.TensorUtils.mapSlices(@min, dim, t);
            val = cell2mat(val);
            idxAlongDim = cell2mat(idxAlongDim);
        end

        function [idxAlongDim, val] = argMax(t, dim)
            % returns the position along dim of the minimum of each
            % position
            [val, idxAlongDim] = Neuropixel.Utils.TensorUtils.mapSlices(@max, dim, t);
            val = cell2mat(val);
            idxAlongDim = cell2mat(idxAlongDim);
        end

        function [coeff, score, latent, tsquared, explained] = pcaMultiDim(t, basisDims, varargin)
            % performs PCA on a matrix where basisDims defines the dimensions over which linear
            % combinations are constructed, and each other dimension
            % represents observations.
            % coeff will be constructed with basisDims reshaped to fill the
            % last dimension, and otherDims as subsequent dimensions

            szOrig = size(t);
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(szOrig, basisDims);
            t = Neuropixel.Utils.TensorUtils.reshapeByConcatenatingDims(t, {otherDims, basisDims});

            % deal gracefully with all nan bases
            validCols = any(~isnan(t), 1);
            [coeff, score, latent, tsquared, explained] = pca(t(:, validCols), varargin{:});
            if ~all(validCols)
                coeff = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(coeff, 1, validCols);
            end
            nC = size(coeff, 2);
            score = reshape(score, [szOrig(otherDims), nC]);
            tsquared = reshape(tsquared, szOrig(otherDims));
        end

        function [coeff, score, latent, tsquared, explained, mu] = pcaAlongDim(t, basisDim, varargin)
            % performs PCA on a matrix where basisDim defines the single dimension over which linear
            % combinations are constructed, and each other dimension
            % represents observations.
            % coeff will be constructed with basisDims reshaped to fill the
            % last dimension, and otherDims as subsequent dimensions

            assert(isscalar(basisDim));
            szOrig = size(t);
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(szOrig, basisDim);
            t = Neuropixel.Utils.TensorUtils.reshapeByConcatenatingDims(t, {otherDims, basisDim});

            % deal gracefully with all nan bases
            validCols = any(~isnan(t), 1);
            [coeff, score, latent, tsquared, explained, mu] = pca(t(:, validCols), varargin{:});

            if ~all(validCols)
                coeff = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(coeff, 1, validCols, 0);
                mu = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(mu', 1, validCols, 0);
            end

            szNew = szOrig;
            nC = size(coeff, 2);
            szNew(basisDim) = nC;
            score = Neuropixel.Utils.TensorUtils.undoReshapeByConcatenatingDims(score, {otherDims, basisDim}, szNew);
            tsquared = reshape(tsquared, szOrig(otherDims));

            mu = Neuropixel.Utils.TensorUtils.orientSliceAlongDims(mu, basisDim);
        end
    end

    methods(Static) % Data manipulation
        function t = centerSlicesOrthogonalToDimension(t, alongDims)
            % takes slices along dimensions specified, computes the mean
            % of the slice and subtracts it. this ensures
            % that the mean of each slice spanning alongDims will be zero
            otherDims = Neuropixel.Utils.TensorUtils.otherDims(size(t), alongDims);
            meanTensor =  Neuropixel.Utils.TensorUtils.nanmeanMultiDim(t, otherDims);
            t = bsxfun(@minus, t, meanTensor);
        end


        function [t, mu] = centerSlicesSpanningDimension(t, alongDims)
            % for each subscript in dimension(s) alongDims, computes the mean
            % along all other dimensions and subtracts it. this ensures
            % that the mean along any slice in alongDims will have zero
            % mean.
            mu =  Neuropixel.Utils.TensorUtils.nanmeanMultiDim(t, alongDims);
            t = bsxfun(@minus, t, mu);
        end

        function t = rescaleIntervalToInterval(t, varargin)
            % Rescale t such that interval `from` is scaled to equal
            % interval `to`.
            %
            % Args:
            %   from (numeric interval) : defaults to [min, max]
            %   to (numeric interval) : defaults to [0 1]
            %

            p = inputParser();
            p.addOptional('from', [nanmin(t(:)), nanmax(t(:))], @isvector);
            p.addOptional('to', [0 1], @isvector);
            p.parse(varargin{:});

            from = p.Results.from;
            to = p.Results.to;

            t = (t - from(1)) / (from(2) - from(1)) * (to(2) - to(1)) + to(1);
        end

        function t = makeNonDecreasing(t, dim)
            if nargin < 2
                dim = Neuropixel.Utils.TensorUtils.firstNonSingletonDim(t);
            end

            % take cumulative max and replace each value that drops below
            % the cumulative max with the cumulative max at that point
            cm = cummax(t, dim);
            mask = t < cm;
            t(mask) = cm(mask);
        end

        function t = makeNonIncreasing(t, dim)
            if nargin < 2
                dim = Neuropixel.Utils.TensorUtils.firstNonSingletonDim(t);
            end

            % take cumulative min and replace each value that rises above
            % the cumulative min with the cumulative min at that point
            cm = cummin(t, dim);
            mask = t > cm;
            t(mask) = cm(mask);
        end

        function reweightedTensor = linearCombinationAlongDimension(t, dim, weightsNewByOld, varargin)
            % looking along dimension dim, reweight the tensor along that
            % dimension by linearly combining. If dim==1 and the tensor
            % is a matrix, this is equivalent to matrix multiplication,
            % reweightedTensor = weightsNewByOld * t.

            p = inputParser();
            p.addParameter('replaceNaNWithZero', false, @islogical); % ignore NaNs by replacing them with zero
            p.addParameter('keepNaNIfAllNaNs', false, @islogical); % when replaceNaNWithZero is true, keep the result as NaN if every entry being combined is NaN
            % on a per-value basis, normalize the conditions by the number of conditions present at that time on the axis
            % this enables nanmean like computations
            p.addParameter('normalizeCoefficientsByNumNonNaN', false, @islogical);
            % requires that the weight matrix be square, the equivalent of adding
            % the identity matrix to the weight matrix, except that this
            % will be added after normalization.
            p.addParameter('addToOriginal', false, @islogical);
            p.parse(varargin{:});

            nOld = size(t, dim);
            assert(size(weightsNewByOld, 2) == nOld, 'Size of weight matrix must have nOld==%d columns', nOld);
            nNew = size(weightsNewByOld, 1);

            if p.Results.addToOriginal
                assert(size(weightsNewByOld, 1) == size(weightsNewByOld, 2), 'Weight matrix must be square for addToOriginal = true');
            end

            sz = size(t);
            newSz = sz;
            newSz(dim) = nNew;

            % put combination dimension first
            pdims = [dim, Neuropixel.Utils.TensorUtils.otherDims(sz, dim)];
            tp = permute(t, pdims);

            % should be nOld x prod(size-t-other-dims)
            tpMat = tp(:, :);

            if p.Results.normalizeCoefficientsByNumNonNaN || p.Results.keepNaNIfAllNaNs
                % count the number of values in each row
                nValidMat = sum(~isnan(tpMat), 1);
            end

            if p.Results.replaceNaNWithZero
                tpMat(isnan(tpMat)) = 0;
            end

            % should be nNew x prod(size-t-other-dims)
            reweightMat = weightsNewByOld * tpMat;

            if p.Results.normalizeCoefficientsByNumNonNaN
                % do the reweighting
                reweightMat = bsxfun(@rdivide, reweightMat, nValidMat);
            end

            if p.Results.addToOriginal
                reweightMat = reweightMat + tpMat;
            end

            if p.Results.keepNaNIfAllNaNs
                % NaN out where nValidMat is zero
                mask = ones(size(nValidMat));
                mask(nValidMat == 0) = NaN;
                reweightMat = bsxfun(@times, reweightMat, mask);
            end

            reweightedTensor = ipermute(reshape(reweightMat, newSz(pdims)), pdims);
        end

        function newTensor = linearCombinationApplyScalarFnAlongDimension(t, dim, logicalNewByOld, fn, varargin)
            % this method is conceptually similar to the linear combination
            % operation, except that it applies fn() to the set of values
            % along dimension dim that are non-zero in logicalNewByOld.

            nOld = size(t, dim);
            assert(size(logicalNewByOld, 2) == nOld, 'Size of weight matrix must have nOld==%d columns', nOld);
            nNew = size(logicalNewByOld, 1);

            % compute min over all trial counts included in each basis
            parts = Neuropixel.Utils.TensorUtils.cellvec(nNew);
            for iNew = 1:nNew
                fnSelect = @(c) fn(c(logicalNewByOld(iNew, :) ~= 0));
                parts{iNew} = Neuropixel.Utils.TensorUtils.mapSlices(fnSelect, dim, t);
            end
            newTensor = cell2mat(cat(dim, parts{:}));  
        end

        function in = assignValueMaskedSelectionAlongDimension(in, dims, mask, value)
            % in = assignValueMaskedSelectionAlongDimension(in, dims, mask, value=NaN)
            %
            % Assign value at each location where mask is true along dimension dim
            % if dims is scalar assign in(..., mask, ...) = value
            % if dims is a vector of dimensions, then mask must be a cell,
            % in which case the assignment is
            % in(..., mask{1}, ..., mask{2}, ...) = value
            %

            if ~exist('value', 'var')
                if iscell(in)
                    value = {[]};
                elseif islogical(in)
                    value = false;
                else
                    value = NaN;
                end
            end
            if iscell(in) && ~iscell(value)
                value = {value};
            end

            assert(isscalar(value), 'Value must be scalar');
            maskByDim = Neuropixel.Utils.TensorUtils.maskByDimCellSelectAlongDimension(size(in), dims, mask);
            in(maskByDim{:}) = value;
        end

        function data = clip(data, clipAt, replace)
            % if replace is empty or not provided, data will be replaced by
            % the correspondign limit.
            % both clipAt can be scalar positive to clip at +/- or [low
            % high] vectors
            if isempty(clipAt)
                return;
            elseif isscalar(clipAt)
                clipAt = [-clipAt clipAt];
            end
            if nargin < 3 || isempty(replace)
                replace = clipAt;
            else
                if isscalar(replace)
                    replace = [-replace replace];
                end
            end

            if ~iscell(data)
                data = clip_I(data);
            else
                for i = 1:numel(data)
                    data{i} = clip_I(data{i});
                end
            end

            function d = clip_I(d)
                d(d < clipAt(1)) = replace(1);
                d(d > clipAt(2)) = replace(2);
            end
        end

    end
end
