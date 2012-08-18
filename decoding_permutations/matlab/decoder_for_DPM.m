% function e = decoder_for_DPM(ftmatrix, tx, q, d, Clist)
% Input:
%     ftmatrix: the output of channel()
%     tx: the q-ary Codeword()._word that is obtained from linear code
%     q:  the alphabet size
%     d:  the minimum distance of the code
%     Clist: the code as a list of lists in the Z_q representation
% Output:
%     The number of symbols in error.
%
function e = decoder_for_DPM(ftmatrix, tx, q, d, Clist)
y = ftmatrix_to_Codeword(ftmatrix);
if q == 2
    s = estimate_binary_DPM(y);
elseif q == 3
    s = estimate_ternary_DPM(y);
else
    error('Not implemented for larger q');
end

num_erasures = length(find(s >= q));
if 2*(hamming_distance(s, tx) - num_erasures) <= d-1 - num_erasures
    e = 0;
    return;
end

% Check if it is closer to some other codeword than to tx
for c = Clist'
    if 2*(hamming_distance(s, c') - num_erasures) <= d-1-num_erasures
        e = hamming_distance(c', tx);
        return;
    end
end

% rx is outside any ball of radius d/2. Decoding failure.
e = hamming_distance(s, tx);
