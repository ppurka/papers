% function txp = binary_to_permutations_by_flip(cw, m)
% Given the codeword cw in Z_{2^m}^n, return the permutation vector
% Example:
%   binary_to_permutations_by_flip([1 1 0 1], 1)
%   ans = 1 2 0 4 3
%
%   binary_to_permutations_by_flip([1 3 2], 2)
%   ans = 0 2 3 4 5 1 6
%
function sigma = binary_to_permutations_by_flip(cw, m)
sigma = 0:m*length(cw);
i = 0;
for c = cw
    j = 0;
    for b = get_binary(c, m)
        if b == 1
            indx = i*m+j+1;
            sigma([indx indx+1]) = sigma([indx+1 indx]); % swap the values
        end
        j = j+1;
    end
    i = i+1;
end

