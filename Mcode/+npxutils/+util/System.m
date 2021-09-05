classdef System
    %SYSTEM System-related utilities
    
    methods (Static)
        function out = userHomeDir
            % Path to this user's home directory.
            %
            % Returns a scalar string containing an absolute path.
            if ispc
                out = getenv('USERPROFILE');
            else
                out = getenv('HOME');
            end
            out = string(out);
        end
    end
    
    methods (Access=private)
        function this = System
            % Constructor is private to prevent instantiation and doc visibility.
        end
    end
end

