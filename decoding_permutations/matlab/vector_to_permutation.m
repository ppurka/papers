% function tx = vector_to_permutation(v, q, Shift)
% map q-ary vector to permutation
% v = the q-ary vector
% q = q of q-ary
% Shift = no. of extra coordinates that will be involved on next symbol
%
% Example:
%   vector_to_permutation([1 2 1], 3)
%   ans = 4   5   3   1   2   6   0
%   vector_to_permutation([1 1 0 1], 2)
%   ans = 3   2   1   4   0
%
function tx = vector_to_permutation(v, q, Shift)
if nargin == 2
    Shift = q-1; % This is ok for q=2,3
end

n = length(v);
tx= generate_permutation(0:(q+Shift*(n-1)-1), v, q, Shift);



function p = generate_permutation(L, u, m, Shift)
if length(u) == 1
    p = qary_shift_left(L, u(1), m);
    return
end
p = generate_permutation(qary_shift_left(L, u(1), m), u(2:end), ...
                        m+Shift, Shift);


% Shift the first m symbols cyclically by b, b in {0,...,q-1}
function p = qary_shift_left(L, b, m)
p = L;
for i = 1:m
    p(i) = mod(L(i)+b, m);
end
