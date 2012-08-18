% function e = decoder_for_binary_DIM(ftmatrix, tx, m, d, Clist)
% Input:
%     ftmatrix: the output of channel()
%     tx: the q-ary Codeword()._word that is obtained from linear code
%     m:  the power in 2^m
%     d:  the minimum distance of the code
%     Clist: the code as a list of lists in the Z_q representation
% Output:
%     The number of symbols in error.
%
% Example:
%   A = eye(8); A(1,1) = 0; z = zeros(1, 7);
%   decoder_for_binary_DIM(A, z, 1, 3, [z; 1 1 1 0 0 0 0])
%   ans = 0
%
%   A(2,1) = 1
%   decoder_for_binary_DIM(A, z, 1, 3, [z; 1 1 1 0 0 0 0])
%   ans = 0
%     
function e = decoder_for_binary_DIM(ftmatrix, tx, m, d, Clist)
e = 0;  % This will hold the number of errors
y = ftmatrix_to_permutation(ftmatrix);
b = estimate_bits(y);
rx= bits_to_symbols(b, m);

num_erasures = length(find(rx >= 2^m));
if 2*(hamming_distance(rx, tx) - num_erasures) <= d-1-num_erasures
    e = 0;
    return;
end

% Check if it is closer to some other codeword than to tx
for c = Clist'   % We must make sure the *columns* are codewords
    if 2*(hamming_distance(rx, c') - num_erasures) <= d-1-num_erasures
        e = hamming_distance(c', tx);
        return
    end
end

% rx is outside any ball of radius d/2. Decoding failure.
e = hamming_distance(rx, tx);
