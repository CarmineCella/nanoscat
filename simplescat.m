% SIMPLESCAT, (c) 2015 carmine e. cella
%
% very simplified implementation of scattering transform
% (dyadic wavelets)

clear all
close all

imunit = i;
downsampling_fac =@(res,j) max(0,floor(j * log2(2)-res-1));


%% load audio

[sig, sr] = audioread ('../../datasets/various_data/5_notes.wav');
sig = sig(1:2^floor(log2(length(sig)))); % truncate to power of 2
sig = sig / norm(sig);
N = length (sig);

%% params
M = 2;
J = 11;

%% compute filters

psi = {};
phi = {};

for res = 0 : floor (log2(N));
    N0 = N / 2^res;
    
    if N0 <= N/2^J;
        break;
    end
    
    for j = 0 : J -1
        v = zeros(1, N0);
        
        sz = N0 /2^j;
        if sz <= N/2^J;
            break;
        end

        s = .5*(1 - cos(2*pi*(0:sz-1)'/(sz)));
        v(1:sz) = s .^1;% hanning-zero
        v(sz) = 0;
        
        psi{res+1}{res+j+1} = v';
 
        if (res+j == J - 1)
            f = zeros(1, N0);
            f(end-(sz/2)+1:end) = v(1:(sz/2));
            f(1:(sz/2))=v((sz/2)+1:sz);
            phi{res+1} = f';
        end
    end
end

%% plot filters


res = 1;
figure
for j = 1:numel(psi{res})
    plot (psi{res}{j});
    hold on
end
plot (phi{res});
%%
accum = zeros (numel(psi{1}{1}), 1);
for i = 1 : numel(psi{1})
    accum = accum + abs (psi{1}{i});
end
accum = accum + abs (phi{1});
plot (accum, 'k')
title ('PSI/PHI  at higher resolution')

%% compute scattering
U = {};
S = {};

U{1}{1}= sig;

for m = 1:M+1
    hindex = 1;
    if m > size(U,2)
        continue;
    end
    
    res = 1;
    for s = 1:numel(U{m})
        sigf = fft (U{m}{s});
        res = (log2(N) - (log2 (length(sigf)))) + 1;
 
        if m<=M
            vindex = 1;
            paths = {};
            for j = s : numel(psi{res})
                ds = 2^(j-s);
                c = abs(ifft(sigf .* psi{res}{j}));
                paths{vindex} = c(1:ds:end);
                vindex = vindex + 1;
            end
            
            for ic=1:numel(paths)
                U{m+1}{hindex}=paths{ic};
                hindex=hindex+1;
            end
        end
        ds = (J - res)^2;
        c = abs (ifft(sigf .* phi{res}));
        if ds > 1
        c = c(1:ds:end);
        end
        S{m}{s} = c; 
    end
end

%% plot S coefficients

for m =1:M+1
    Out = S{m};
    
    maxlen = length (Out{1});
    nlambdas = numel(Out);
    
    Scoeff = zeros(nlambdas, maxlen);
    
    for i = 1:nlambdas        
        %Scoeff(i, 1:length(Out{i})) = Out{i};
        itp = interpft(Out{i}, maxlen);
        Scoeff(i, :) = itp;
    end
    
    figure
    imagesc (Scoeff)
    title (sprintf ('Order %d', m - 1))
end

% eof

