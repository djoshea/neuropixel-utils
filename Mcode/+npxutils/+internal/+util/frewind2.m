function frewind2(fid)
% A version of frewind that errors on failure
npxutils.internal.util.fseek(fid, 0, 'bof');
end