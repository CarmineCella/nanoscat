function [sig, N, len] = nanoscat_load (filename)
[y, sr] = audioread (filename);
len=length(y);
N = 2^(floor(log2(len)) + 1);
sig = zeros(N, 1);
sig(1:len) = y;
end
