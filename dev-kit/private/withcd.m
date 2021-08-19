function out = withcd(dir)
% Temporarily change to a new directory
origDir = pwd;
cd(dir);
out.RAII = onCleanup(@() cd(origDir));
end