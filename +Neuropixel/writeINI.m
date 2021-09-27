function writeINI(fname, data, varargin)
    p = inputParser();
    p.addParameter('prefixTildeFields', strings(0, 1), @isstring); % these fields should have a tilde before them, ignoring this broke certain Python tools like Neo
    p.addParameter('quoteStrings', false, @islogical);
    p.parse(varargin{:});
    prefixTildeFields = string(p.Results.prefixTildeFields);
    quoteStrings = p.Results.quoteStrings;
    
    fid = fopen(fname,'w');
    if fid == -1
        error('Could not open %s for writing', fname);
    end
    
    flds = fieldnames(data);
    
    for iF = 1:numel(flds)
        valStr = convertVal(data.(flds{iF}), quoteStrings);
        
        if ismember(flds{iF}, prefixTildeFields)
            fld = "~" + string(flds{iF});
        else
            fld = flds{iF};
        end
        fprintf(fid, '%s=%s\n', fld, valStr);
    end
    
    fclose(fid);
end

function valStr = convertVal(val, quoteStrings)
    if ischar(val)
        if isempty(val)
            valStr = '';
        elseif quoteStrings
            valStr = ['''' strtrim(val) ''''];
        else
            valStr = strtrim(val);
        end
        
    elseif isstring(val)
        if isscalar(val)
            if quoteStrings
                valStr = "'" + strtrim(val) + "'";
            else
                valStr = strtrim(val);
            end
        else
            valStrCell = cellstr(val);
            valStr = ['{' strjoin(valStrCell, ',') '}'];
        end
    elseif iscell(val)
        valStrCell = cellfun(@convertVal, val(:), 'UniformOutput', false);
        valStr = ['{' strjoin(valStrCell, ',') '}'];
        
        
    elseif islogical(val)
        if val
            valStr = 'true';
        else
            valStr = 'false';
        end
        
    elseif isscalar(val)
            if isinteger(val) || abs(round(val) - val) < val / 10^9
                valStr = sprintf('%d', round(val));
            else
                valStr = sprintf('%g', val);
            end
            
    
    elseif isvector(val)
        valStrCell = arrayfun(@convertVal, val(:), 'UniformOutput', false);
        valStr = strjoin(valStrCell, ',');
        
    else
        valStr = mat2str(val);
    end
end