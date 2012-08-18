% function ftmatrix = generate_ftmatrix(word, q)
% Generate the frequency-time matrix of the codeword "word"
% Inputs:
%   word:   the codeword
%   q:      the alphabet size
%
% Output:
%   ftmatrix: the frequency-time matrix
%
% Example:
%   generate_ftmatrix([0 1 2], 2)
%   ans =    1   0   0
%            0   1   0
%
function ftmatrix = generate_ftmatrix(word, q)
ftmatrix = zeros(q, length(word));
i = 1;
for w = word
    if w < q
        ftmatrix(w+1, i) = 1;
    end
    i = i+1;
end
