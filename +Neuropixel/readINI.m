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
    
    vald = str2double(val);
    if strcmpi(val, 'NaN')
        val = NaN;
    elseif ~any(ismember({','}, val)) && ~isnan(vald)
        val = vald;
    elseif strcmpi(val, 'false')
        val = false;
    elseif strcmpi(val, 'true')
        val = true;
    elseif val(1) == '''' && val(end) == ''''
        val = val(2:end-1);
    end
end