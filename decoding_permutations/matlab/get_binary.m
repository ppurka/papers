% get the binary representation of the symbol c as a list of length m.
% Example:
%   get_binary(2, 2)
%   ans = 1 0
function binary = get_binary(c, m)
binary = zeros(1, m);
for i = m:-1:1
    binary(i) = mod(c, 2);
    c = (c - binary(i))/2;
end
