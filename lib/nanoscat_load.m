function [sig, N] = nanoscat_load (filename)
[y, sr] = audioread (filename);
N = 2^(floor(log2(length(y))) + 1);
sig = zeros(N, 1);
sig(1:length(y)) = y;
end
