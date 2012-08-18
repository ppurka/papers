% function cw = ftmatrix_to_permutation(ftmatrix)
% This function converts a frequency-time matrix to a Codeword. At time
% instance i, it looks only on the frequencies which are less than i+1.
% The code is identical to the other function ftmatrix_to_Codeword except
% for the step where we check for an erasure (the last for loop).
% The steps followed are:
% 1. If there is a list of all 1's in a particular frequency, then it is
%    set to an all 0 list
% 2. If there is a list of all 1's in a particular time, then it is set
%    to an all 0 list
% 3. If there is all 0 time coordinate then it is set to q (an erasure
%    symbol)
% 4. If there are more than two 1's at time instance i then if one of the
%    1's is at (i+1,i)-th frequency-time position, then the coordinate is
%    assumed to be i+1. If there are more than two 1's in (j,i)-th
%    frequency-time positions with j<i then the first frequency index of
%    the 1 is assumed as the coordinate.
% 5. Otherwise it returns the index of the 1 in a particular time.
%
% Output is a Codeword with symbols from {0,...,q} with q denoting an
% erasure symbol.
%
% Example:
%   A = eye(7); A(2,:) = 1; A(:,[3, 4]) = 1; A(3,6) = 1;
%   ftmatrix_to_permutation(A)
%   ans = 0   7   7   7   4   2   6
%
%   A = generate_ftmatrix([1   2   0   3   5   4   7   6], 8);
%   ftmatrix_to_permutation(A)
%   ans = 1   2   0   3   5   4   7   6
%
function cw = ftmatrix_to_permutation(ftmatrix)
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
for f = ftmatrix        % iterate over the columns of the q x n matrix
    E = find(f == 1);
    %disp(["DEBUG: E after find: ", num2str(E')])
    E = E(E <= i+1);    % have to be careful with indices. crazy matlab.
    %disp(["DEBUG: E after <=i+1: ", num2str(E')])
    flag = length(find(E==i+1));    % This is set to 1 if i is in E
    %disp(["DEBUG: flag = ", num2str(flag)])
    if flag
        cw(i) = i;
    elseif length(E) >= 1
        cw(i) = E(1)-1;
    else
        cw(i) = q;
    end
    i = i+1;
end
