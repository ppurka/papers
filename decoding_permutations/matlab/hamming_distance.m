% function d = hamming_distance(x, y)
% Compute the Hamming distance between two vectors x and y
%
function d = hamming_distance(x, y)

if length(x) ~= length(y)
    error('The lengths of the input vectors must be same');
end

d = sum(x ~= y);