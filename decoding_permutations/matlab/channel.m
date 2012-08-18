% function ftmatrix = channel(word, q, impulse_noise, nb_noise, prob)
% input:
%     word: a Codeword. This is the transmitted word.
%
%     q: the alphabet size
% 
%     impulse_noise: impulse noise. If the parameter is a positive integer
%     then it acts as a worst-case channel. If the parameter is between
%     0 and 1 then it acts as a random channel.
% 
%     nb_noise: narrow band noise. If the parameter is a positive integer
%     then it acts as a worst-case channel. If the parameter is between
%     0 and 1 then it acts as a random channel.
% 
%     prob: probability with which the additive white noise will be
%     generated at any given time instance. This must be strictly in the
%     interval [0,1]. The additive noise is generated for each frequency
%     and time and so it will insert a 1 in the frequency-time matrix
%     according to the frequency and time it corresponds to.
% 
% output: a freq-time matrix containing codeword + impulse noise
% + narrow-band noise + random noise
% 
% In the frequency-time matrix, the alphabet is {0,1}. With 0 being
% reserved for the case when no errors occur. The position of 1 gives the
% frequency along the row and the time along the column. For instance,
% the matrix [0 0;1 0;0 1] corresponds to the codeword (1,2).
%
% Example:
%   word = 0:100;
%   sum(sum(channel(word, length(word), 0, 0, 0.1) ~= ...
%       generate_ftmatrix(word, length(word))))/(length(word)^2)
%   ans =  0.10401 % random answer
%
%   sum(sum(channel(word, length(word), 0, 0, 0.2) ~= ...
%       generate_ftmatrix(word, length(word))))/(length(word)^2)
%   ans =  0.1992 % random answer
%
function ftmatrix = channel(word, q, impulse_noise, nb_noise, prob)
if nargin < 5, prob = 0;            end
if nargin < 4, nb_noise = 0;        end
if nargin < 3, impulse_noise = 0;   end

n = length(word);

% Add random noise. bsc simulates the binary symmetric channel. It flips
% every bit with probability prob.
ftmatrix = bsc(generate_ftmatrix(word, n), prob);

% TODO: Add code for impulse noise and nb_noise.
%       right now, we are using only background noise, so we don't need
%       this

% OCTAVE SPECIFIC CODE. Comment this function when using MATLAB
function ftmatrix = bsc(ftmatrix, prob)
M = (rand(size(ftmatrix)) <= prob);
ftmatrix = mod(ftmatrix + M, 2);
