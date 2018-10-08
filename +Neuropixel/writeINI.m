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
        valStr = ['''' strtrim(val) ''''];
    elseif islogical(val)
        if val
            valStr = 'true';
        else
            valStr = 'false';
        end
    elseif isscalar(val)
            if isinteger(val)
                valStr = sprintf('%d', val);
            else
                valStr = sprintf('%g', val);
            end
    else
        valStr = [mat2str(val)];
    end
end