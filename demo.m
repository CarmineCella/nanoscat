% nanoscat demo, (c) 2015 carmine e. cella
%
% very simplified implementation of scattering transform for 1D signals
% for pedagogical purpouses (dyadic wavelets)
%
clear all
close all

addpath ('lib');

%% params
M = 2; % orders
J = 12; % maximal scale

%% load and zero pad audio
[sig, N] = nanoscat_load ('samples/drum1_90.wav');
sig = sig / norm(sig); % normalization

assert (J <= log2(N));

%% compute filters
[psi, phi, lp] = nanoscat_filters (N, J);

%% plot filters
nanoscat_display_filters (psi, phi, lp);

%% compute scattering
[S, U] = nanoscat_compute (sig, psi, phi, M);

%% format and plot S coefficients
scat = nanoscat_format (S, [1:M+1]); % creates a matrix with all coefficients

figure
imagesc (scat);
title ('Scattering coefficients (all orders)');

% eof

