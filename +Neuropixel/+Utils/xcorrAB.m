function C = xcorrAB(A, B, maxlag)
% like xcorr but computed between all the columns of A and the columns of B.
% A and B must have the same number of rows
%
% C is 2*maxlag-1 x size(A, 2)*size(B, 2). 
% You may want to call reshape(C, [], size(A, 2), size(B, 2)) to make the 
% output more straightforward to use

    [m,n] = size(A);
    assert(size(B, 1) == m, 'A and B must have same number of rows');
    
    if nargin < 3
        maxlag = m - 1;
    end
    mxl = min(maxlag,m - 1);

    m2 = findTransformLength(m);

    XA = fft(A,m2,1);
    XB = fft(B,m2,1);
    C = reshape(XA,m2,1,n).*conj(XB(:,:));

    % Call IFFT and force real output if x is real.
    c1 = ifft(C,[],1,'symmetric');
    % c1 is M2-by-N-by-N.
    % Keep only the lags we want, and move the negative lags before the
    % positive lags. Also flatten the result to 2-D.
    C = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];

end

function m = findTransformLength(m)
    m = 2*m;
    while true
        r = m;
        for p = [2 3 5 7]
            while (r > 1) && (mod(r, p) == 0)
                r = r / p;
            end
        end
        if r == 1
            break;
        end
        m = m + 1;
    end
end