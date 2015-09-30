% SIMPLESCAT, (c) 2015 carmine e. cella
%
% very simplified implementation of scattering transform
% (dyadic wavelets)

clear all
close all

imunit = i;

%% load audio

[sig, sr] = audioread ('../../datasets/various_data/5_notes.wav');
sig = sig(1:2^floor(log2(length(sig)))); % truncate to power of 2
sig = sig / norm(sig);
N = length (sig);

%% params
M = 2;
J = 11;

%% compute filters
res = 0;
psi = {};
phi = {};

for res = 0 : floor (log2(N));
    N0 = N / 2^res;
    
    for j = 0 : J - 1
        v = zeros(1, N0);
        
        sz = N0/2^j;
        if sz <= 2
            continue
        end
        v(1:sz) = .5*(1 - cos(2*pi*(0:sz-1)'/(sz))); % hanning-zero
        v(sz) = 0;
        
        psi{res+1}{res+j+1} = v';
        if j == 0
            phi{res+1} = 1 ./ (1 + v');
        end
    end
end

%% plot filters
%for res = 1 : numel(psi)
    res = 1;
    figure
    for j = 1:numel(psi{res})
        plot (psi{res}{j});
        hold on
    end
%end
title ('PSI (higher resolution)')

figure
for res = 1:numel(phi)
    plot (phi{res});
    hold on
end
title ('PHI')

%% compute scattering
U = {};
S = {};

U{1}{1}= sig;

for m = 1:M+1
    hindex = 1;
    if m > size(U,2)
        continue;
    end
    
    for s = 1:numel(U{m})
        sigf = fft (U{m}{s});
        
        if m<=M
            vindex = 1;
            paths = {};
            for j = s: numel(psi{1})-1
                paths{vindex} = abs(ifft(sigf .* psi{1}{j})); % FIXME: downsampling is missing5
                vindex = vindex + 1;
            end
            
            for ic=1:numel(paths)
                U{m+1}{hindex}=paths{ic};
                hindex=hindex+1;
            end
        end
        
        S{m}{s} = abs (ifft(sigf .* phi{1}));
    end
    
end

%% plot S coefficients

for m =1:M+1
    Out = S{m};
    
    maxlen = length (Out{1});
    nlambdas = numel(Out);
    
    Scoeff = zeros(nlambdas, maxlen);
    
    for i = 1:nlambdas        
        Scoeff(i, 1:length(sig)) = Out{i};
        %itp = interpft(Out{i}.signal, maxlen);
        %Scoeff(i, :) = itp;
    end
    
    figure
    imagesc (Scoeff)
    title (sprintf ('Order %d', m - 1))
end

% eof

