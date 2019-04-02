function success = makeSymLink(src, link, relative)
% makeSymLink(src, linkDest)

    if nargin < 3
        relative = false;
    end

    src = resolveSymLink(Neuropixel.Utils.GetFullPath(src));
    if ~exist(src, 'file')
        warning('Source is a dangling symlink, not creating');
    end
    
    link = Neuropixel.Utils.GetFullPath(link);
    Neuropixel.Utils.mkdirRecursive(fileparts(link));
    % cant use exist on symlinks since it tests for existence of the file it points to
    sWarn = warning('off', 'MATLAB:DELETE:FileNotFound');
    try 
        delete(link);
    catch
    end
    warning(sWarn);
    
    if relative
        src = Neuropixel.Utils.relativepath(src, fileparts(link));
    end
    
    cmd = sprintf('ln -s "%s" "%s"', src, link);
    [status, output] = unix(cmd);
    
    if status
        fprintf('Error creating symlink: \n');
        fprintf('%s\n', output);
    end

    success = ~status;
end

function target = resolveSymLink(file)

    if ismac
        target = file;
        return;
    else
        target = getResolved(file);
    end
    
    return;

    if ~exist(file, 'file')
        % try recursively on its parent
        [parent leaf ext] = fileparts(file);
        parent = resolveSymLink(parent);
        target = fullfile(parent, [leaf ext]);
    else
        target = getResolved(file);
    end

    return;
    
    function result = getResolved(file)
        cmd = sprintf('readlink -m "%s"', Neuropixel.Utils.GetFullPath(file));
        [status, result] = system(cmd);
        if status || isempty(result)
            fprintf(result);
            error('Error resolving sym link');
        end 
        
        NEWLINE = 10;
        if result(end) == NEWLINE
            result = result(1:end-1);
        end
    end
end
