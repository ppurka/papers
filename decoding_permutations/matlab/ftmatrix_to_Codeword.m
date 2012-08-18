% function cw = ftmatrix_to_Codeword(ftmatrix)
% This function converts a frequency-time matrix to a Codeword.
% The steps followed are:
% 1. If there is a list of all 1's in a particular frequency, then it is
%    set to an all 0 list
% 2. If there is a list of all 1's in a particular time, then it is set
%    to an all 0 list
% 3. If there is all 0 all time coordinate or more than one 1 in a time
%    coordinate then it is set to q (an erasure symbol)
% 4. Otherwise it returns the index of the 1 in a particular time.
%
% Output is a Codeword with symbols from {0,...,q} with q denoting an
% erasure symbol.
%
% Example:
%   ftmatrix_to_Codeword(generate_ftmatrix(0:2, 2))
%   ans = 0   1   2
%   A = eye(7); A(2,:) = 1; A(:,[3, 4]) = 1; A(3,6) = 1;
%   ftmatrix_to_Codeword(A)
%   ans = 0   7   7   7   4   7   6
%
function cw = ftmatrix_to_Codeword(ftmatrix)
[q n] = size(ftmatrix); % q = # rows, n = # cols

% find all the symbols with narrow band errors
narrow_band_err = zeros(1, q);
o = ones(n, 1);
i = 1;
for l = ftmatrix' % iterate over the columns of the n x q matrix
    if l == o
        narrow_band_err(i) = 1;
    end
    i = i+1;
end
narrow_band_err = find(narrow_band_err == 1);

% find all the impulse noise errors
impulse_err = zeros(1, n);
o = ones(q, 1);
i = 1;
for l = ftmatrix
    if l == o
        impulse_err(i) = 1;
    end
    i = i+1;
end
impulse_err = find(impulse_err == 1);

% Delete the narrow band and impulse errors
ftmatrix(narrow_band_err, :) = 0;
ftmatrix(:, impulse_err)     = 0;

% In this code ftmatrix is not transposed yet
cw = zeros(1, n);
i  = 1;
for f = ftmatrix
    e = find(f == 1);
    if length(e) == 1
        cw(i) = e(1)-1;
    else
        cw(i) = q;
    end
    i = i+1;
end
