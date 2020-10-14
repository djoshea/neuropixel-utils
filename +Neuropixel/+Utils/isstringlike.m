function tf = isstringlike(str)
    tf = ischar(str) || isstring(str) || iscellstr(str);
end