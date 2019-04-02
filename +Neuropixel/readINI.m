function out = readINI(fname)
    out = struct();
    
    if ~exist(fname, 'file')
        error('Configuration file %s not found', fname);
    end
    f = fopen(fname,'r');
    while ~feof(f)
        s = strtrim(fgetl(f));
        if isempty(s) || s(1)==';' || s(1)=='#'
            continue;
        end
        
        [par, val] = parseLine(s);
        out.(par) = val;
    end
    
    fclose(f);
end

function [par, val] = parseLine(s)
    [par,val] = strtok(s, '=');

    par = strtrim(par);
    if par(1) == '~'
        par(1) = [];
    end
    par = strtrim(par);
    par = genvarname(par);
    
    val = strtrim(val);
    if val(1) == '='
        val(1) = [];
    end
    val = strtrim(val);
    if isempty(val)
        return;
    end
    
    [vald, valid] = str2vector(val);
    
    if all(valid)
        val =vald;
    else
        if strcmpi(val, 'NaN')
            val = NaN;
        elseif ~any(ismember({','}, val)) && ~any(isnan(vald))
            val = vald;
        elseif strcmpi(val, 'false')
            val = false;
        elseif strcmpi(val, 'true')
            val = true;
        elseif val(1) == '''' && val(end) == ''''
            val = val(2:end-1);
        end
    end
end

function [vec, valid] = str2vector(str, varargin)
% parses a vector specified as a string that looks like
% '[1,2,3,4]', '[1 2 3 4]', '1 2 3 4', or '1,2,3,4'
% 
%  vec is the parsed vector or [] if valid == false
%  if the string is blank or the vector doesn't contain any values
%  vec will == emptyValue, which is by default NaN (not []) and can be overriden
%  using str2vector(str, 'emptyValue', emptyValue);
%
%  if str is a cell array of strings, vec will be a cell array of vectors and 
%  valid will be a logical array of the same size

emptyValue = NaN;
invalidValue = NaN;
str = strtrim(str);

if isempty(str)
    vec = emptyValue;
    valid = true;
    return;
end

% figure out which delimiter we're using to separate the values
% semicolon, comma, or space
if any(str == ';')
    delim = ';';
    multDelimsAsOne = false;
    orientFn = @Neuropixel.Utils.makecol;
elseif any(str == ',')
    delim = ',';
    multDelimsAsOne = false;
    orientFn = @Neuropixel.Utils.makerow;
else
    delim = ' ';
    multDelimsAsOne = true;
    orientFn = @Neuropixel.Utils.makerow;
end

asCell = false;
% strip surrounding [] and ''
str = strip(str, '''');
str = strip(str, '[');
str = strip(str, ']');
if ~isempty(str) && str(1) == '{'
    asCell = true;
    str = str(2:end);
end
if ~isempty(str) && str(end) == '}'
    str = str(1:end-1);
end
str = strtrim(str);

if isempty(str)
    vec = emptyValue;
    valid = true;
    return;
end

results = textscan(str, '%s', 'Delimiter', delim, 'MultipleDelimsAsOne', multDelimsAsOne);
tokens = results{1};

ignoreMask = cellfun(@(x) isempty(x) || strcmp(x, 'NaN'), tokens);

if ~asCell
    vec = str2double(tokens);
    valid = ~any(isnan(vec) & ~ignoreMask);
    vec = orientFn(vec);
else
    vec = cell(1, numel(tokens));
    for iT = 1:numel(tokens)
        test = str2double(tokens{iT});
        if isnan(test) && ~strcmp(tokens{iT}, 'NaN')
            % treat as string, strip quotes
            vec{iT} = strip(strip(tokens{iT}, ''''), '"');
        else
            vec{iT} = test;
        end
    end
            
    valid = true;
end


end