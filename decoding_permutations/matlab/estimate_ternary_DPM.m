% function s = estimate_ternary_DPM(y)
% Estimate the ternary vector given the received permutation vector y
%
% Example:
%   estimate_ternary_DPM([4 5 3 1 2 6 0])
%   ans = 1   2   1
%
function s = estimate_ternary_DPM(y)
e = 3;          % the erasure symbol for ternary
N = length(y);  % the alphabet size of permutation space
L =  ones(1,N)*(-1); % the list of non-erased permutation symbols
n = (N-1)/2;    % the length of ternary vector
s = zeros(1,n); % this will hold the estimated ternary symbols

for j = 1:n
    J = 2*(j-1);
    if J-1 >= 0 && y(J) < N
        L(J) = y(J);
    end
    if y(J+1) < N
        L(J+1) = y(J+1);
    end

    J = 2*j;
    if y(J+1) >= N || y(J) >= N
        s(j) = e;
    else
        t = zeros(1, 3);
        for yl = L(L>-1)
            p = [yl-y(J) yl-y(J+1)];
            if      p(1) < 0 && p(2) < 0, t(1) = t(1)+1;
            elseif  p(1) > 0 && p(2) > 0, t(3) = t(3)+1;
            elseif  p(1) < 0 && p(2) > 0, t(2) = t(2)+1;
            end
        end

        if t == zeros(1, 3)
            s(j) = e;
        else
            [~, argmaxtmp] = max(t);
            s(j) = argmaxtmp-1;
        end
    end
end
