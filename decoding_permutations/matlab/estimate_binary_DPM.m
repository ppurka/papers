% function b = estimate_binary_DPM(y)
% Estimate the binary vector given the received permutation vector y
%
% Example:
%   estimate_binary_DPM(vector_to_permutation([1 1 0 1], 2))
%   ans = 1   1   0   1
%
function b = estimate_binary_DPM(y)
b = zeros(1, length(y)-1);
e = 2;          % The erasure bit
L = ones(1,length(y))*(-1); % Unlike the algorithm in paper, these are not the
                % indices but the actual symbols.
N = length(y);  % This is the symbol size of permutation space
i = 1;
for r = y(2:end)
    % Note that r points to y[i+1] and erasure => y[i] >= N
    if y(i) < N % && length(L) < l0 % This can be used to restrict |L|
        L(i) = y(i);
    end
    if r >= N   % the permutations are from 0:N-1. So this denotes erasure
        b(i) = e;
    else
        t = 0;
        for yl = L(L>-1)
            if yl - r > 0
                t = t+1;
            else
                t = t-1;
            end
        end
        if t > 0
            b(i) = 1;
        elseif t < 0
            b(i) = 0;
        else
            b(i) = e;
        end
    end
    i = i+1;
end
