% function b = estimate_bits(y)
% Given the received vector y, estimate the bits. This does not convert to
% the 2^m ary symbols.
%
% Example:
%   estimate_bits([1 2 0 4 3])
%   ans = 1   1   0   1
%   A = eye(7); % identity permutation, so we send all zero in bits
%   A(2,:) = 1; A(:,[3, 4]) = 1; A(3,6) = 1; % the errors
%   ftmatrix_to_permutation(A)
%   ans = 0   7   7   7   4   2   6
%   estimate_bits([ 0   7   7   7   4   2   6])
%   ans = 0   2   2   2   0   0
function b = estimate_bits(y)
N = length(y);
b = zeros(1, N);
i = 1;
for yi = y
    if yi == i
        b(i) = 1;
    elseif yi < i
        b(i) = 0;
    else
        b(i) = 2;
    end
    i = i+1;
end
b = b(1:N-1); % The last symbol is not needed
