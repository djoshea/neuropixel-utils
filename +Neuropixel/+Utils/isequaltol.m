function tf = isequaltol(a,b,tol)
    if isnumeric(a) && isnumeric(b)
        if ~isequal(size(a), size(b))
            tf = false;
            return;
        end
        if nargin < 3
            tol = 1e-6;
        end
        tf = all(abs(a(:) - b(:)) <= tol);
    else
        tf = isequal(a, b);
    end 
end