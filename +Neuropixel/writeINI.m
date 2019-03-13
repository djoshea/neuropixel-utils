function writeINI(fname, data)
    fid = fopen(fname,'w');
    if fid == -1
        error('Could not open %s for writing', fname);
    end
    
    flds = fieldnames(data);
    
    for iF = 1:numel(flds)
        valStr = convertVal(data.(flds{iF}));
        fprintf(fid, '%s=%s\n', flds{iF}, valStr);
    end
    
    fclose(fid);
end

function valStr = convertVal(val)
    if ischar(val)
        if isempty(val)
            valStr = '';
        else
            valStr = ['''' strtrim(val) ''''];
        end
        
    elseif isstring(val)
        valStrCell = cellstr(val);
        valStr = ['{' strjoin(valStrCell, ',') '}'];
        
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