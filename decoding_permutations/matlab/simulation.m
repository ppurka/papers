%-------------------------------------------------------%
%                                                       %
%       Run this block before every simulation          %
%                                                       %
%-------------------------------------------------------%

OCTAVE          = 0;
o               = ones(1, 9);
num_tx          = [o*10^3 o*10^4 o*10^5];
prob_list       = [ 9*10^(-1):-10^(-1):10^(-1) ...
                    9*10^(-2):-10^(-2):10^(-2) ...
                    9*10^(-3):-10^(-3):10^(-3) ];
%num_tx          = [o];
%prob_list       = [ 9*10^(-3):-10^(-3):10^(-3) ];



%-------------------------------------------------------%
%                                                       %
%           Simulation for DIM over GF(2)               %
%                                                       %
%-------------------------------------------------------%
% Set the random seed
if OCTAVE
    rand('seed', 0);                    % OCTAVE code
else
    reset(RandStream.getDefaultStream); % MATLAB code
end

num_symbol_err  = zeros(1, length(prob_list));
[Clist d]       = load_BCH_code('DIM_2'); % every row is a codeword
[M n]           = size(Clist);
m = 1;  % This is the exponent in 2^1
q = 2^m;
i = 1;
for prob = prob_list
    for j = 1:num_tx(i)
        tx  = Clist(randi(M), :);
        txp = binary_to_permutations_by_flip(tx, m);

        rx  = channel(txp, q, 0, 0, prob);
        num_symbol_err(i) = num_symbol_err(i) + ...
                            decoder_for_binary_DIM(rx, tx, m, d, Clist);
    end
    num_symbol_err(i) = num_symbol_err(i)/(num_tx(i)*n);
    fprintf(1, '(prob, sym err rate) = (%f, %f)\n', prob, num_symbol_err(i));
    if OCTAVE, fflush(stdout); end                  % Enable in OCTAVE
    i = i+1;
end
SER_2_DIM = num_symbol_err;
disp(['SER_2_DIM = ', num2str(SER_2_DIM)]);
if OCTAVE, fflush(stdout); end                      % Enable in OCTAVE
figure();
loglog(prob_list, SER_2_DIM, '*');
title('SER 2 DIM - map \Pi_0');

%-------------------------------------------------------%
%                                                       %
%           Simulation for DIM over GF(4)               %
%                                                       %
%-------------------------------------------------------%
% Set the random seed
if OCTAVE
    rand('seed', 0);                    % OCTAVE code
else
    reset(RandStream.getDefaultStream); % MATLAB code
end

num_symbol_err  = zeros(1, length(prob_list));
[Clist d]       = load_BCH_code('DIM_4'); % every row is a codeword
[M n]           = size(Clist);
m = 2;  % This is the exponent in 2^2
q = 2^m;
i = 1;
for prob = prob_list
    for j = 1:num_tx(i)
        tx  = Clist(randi(M), :);
        txp = binary_to_permutations_by_flip(tx, m);

        rx  = channel(txp, q, 0, 0, prob);
        num_symbol_err(i) = num_symbol_err(i) + ...
                            decoder_for_binary_DIM(rx, tx, m, d, Clist);
    end
    num_symbol_err(i) = num_symbol_err(i)/(num_tx(i)*n);
    fprintf(1, '(prob, sym err rate) = (%f, %f)\n', prob, num_symbol_err(i));
    if OCTAVE, fflush(stdout); end                  % Enable in OCTAVE
    i = i+1;
end
SER_4_DIM = num_symbol_err;
disp(['SER_4_DIM = ', num2str(SER_4_DIM)]);
if OCTAVE, fflush(stdout); end                      % Enable in OCTAVE
figure();
loglog(prob_list, SER_4_DIM, '*');
title('SER 4 DIM - map  \Pi_1 ');


%-------------------------------------------------------%
%                                                       %
%           Simulation for DPM over GF(2)               %
%                                                       %
%-------------------------------------------------------%
% Set the random seed
if OCTAVE
    rand('seed', 0);                    % OCTAVE code
else
    reset(RandStream.getDefaultStream); % MATLAB code
end

num_symbol_err  = zeros(1, length(prob_list));
[Clist d]       = load_BCH_code('DPM_2'); % every row is a codeword
[M n]           = size(Clist);
q = 2;
i = 1;
for prob = prob_list
    for j = 1:num_tx(i)
        tx  = Clist(randi(M), :);
        txp = vector_to_permutation(tx, q);

        rx  = channel(txp, q, 0, 0, prob);
        num_symbol_err(i) = num_symbol_err(i) + ...
                            decoder_for_DPM(rx, tx, q, d, Clist);
    end
    num_symbol_err(i) = num_symbol_err(i)/(num_tx(i)*n);
    fprintf(1, '(prob, sym err rate) = (%f, %f)\n', prob, num_symbol_err(i));
    if OCTAVE, fflush(stdout); end                  % Enable in OCTAVE
    i = i+1;
end
SER_2_DPM = num_symbol_err;
disp(['SER_2_DPM = ', num2str(SER_2_DPM)]);
if OCTAVE, fflush(stdout); end                      % Enable in OCTAVE
figure();
loglog(prob_list, SER_2_DPM, '*');
title('SER 2 DPM - map  \Pi_2 ');

%-------------------------------------------------------%
%                                                       %
%           Simulation for DPM over GF(3)               %
%                                                       %
%-------------------------------------------------------%
% Set the random seed
if OCTAVE
    rand('seed', 0);                    % OCTAVE code
else
    reset(RandStream.getDefaultStream); % MATLAB code
end

num_symbol_err  = zeros(1, length(prob_list));
[Clist d]       = load_BCH_code('DPM_3'); % every row is a codeword
[M n]           = size(Clist);
q = 3;
i = 1;
for prob = prob_list
    for j = 1:num_tx(i)
        tx  = Clist(randi(M), :);
        txp = vector_to_permutation(tx, q);

        rx  = channel(txp, q, 0, 0, prob);
        num_symbol_err(i) = num_symbol_err(i) + ...
                            decoder_for_DPM(rx, tx, q, d, Clist);
    end
    num_symbol_err(i) = num_symbol_err(i)/(num_tx(i)*n);
    fprintf(1, '(prob, sym err rate) = (%f, %f)\n', prob, num_symbol_err(i));
    if OCTAVE, fflush(stdout); end                  % Enable in OCTAVE
    i = i+1;
end
SER_3_DPM = num_symbol_err;
disp(['SER_3_DPM = ', num2str(SER_3_DPM)]);
if OCTAVE, fflush(stdout); end                      % Enable in OCTAVE
figure();
loglog(prob_list, SER_3_DPM, '*');
title('SER 3 DPM - map  \Pi_3 ');
