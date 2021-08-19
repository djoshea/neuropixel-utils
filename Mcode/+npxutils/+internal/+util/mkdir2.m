function mkdir2(varargin)
% A version of mkdir that raises error on failure

[ok,msg] = mkdir(varargin{:});
if ~ok
  if nargin == 1
    target = varargin{1};
  else
    target = fullfile(varargin{:});
  end
  error('Failed creating directory "%s": %s', target, msg);
end

end
