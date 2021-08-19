function out = size2str(sz)
%SIZE2STR Format an array size for display
%
% out = size2str(sz)
%
% Sz is an array of dimension sizes, in the format returned by SIZE.
%
% Examples:
%
% size2str(size(magic(3)))

strs = repmat(string(missing), size(sz));
for i = 1:numel(sz)
	strs(i) = sprintf("%d", sz(i));
end

out = strjoin(strs, '-by-');

end
