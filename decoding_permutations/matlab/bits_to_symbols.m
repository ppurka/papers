% function c = bits_to_symbols(b, m)
% Convert a binary list of length nm into a 2^m-ary list of length n by
% combining every m successive bits. The layout of the bits is MSB first.
% If there is an "erasure" (value 2) among the m consecutive bits, then
% simply make the symbol an erasure.
%
function c = bits_to_symbols(b, m)
if m == 1
    c = b;
    return;
end

% now m > 1
c = zeros(1, floor(length(b)/m));
i = 0;
j = 0;
q = 2^m;
s = 0;
twolist = 2.^(m-1:-1:0);
for bj = b
    if bj > 1
        s = q;
    else
        s = s+twolist(j+1)*bj;
    end
    j = j+1;
    if mod(j, m) == 0
        c(i+1) = s; % If s is an erasure, this value is >= q
        i = i+1;
        j = 0;
        s = 0;
    end
end
