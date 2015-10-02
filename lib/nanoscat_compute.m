function [S, U] = nanoscat_compute (sig, psi, phi, M)
U = {};
S = {};

U{1}{1}= sig;
log2N = log2 (length (psi{1}{1})); % maximal length
J = numel(psi);

for m = 1:M+1
    hindex = 1;
    if m > size(U,2)
        continue;
    end
    
    res = 1;
    for s = 1:numel(U{m})
        sigf = fft (U{m}{s});
        res = (log2N - (log2 (length(sigf)))) + 1;
        
        if m<=M
            vindex = 1;
            paths = {};
            for j = s : numel(psi{res})
                ds = 2^(j-s);
                c = abs(ifft(sigf .* psi{res}{j}));
                if res > 2
                    paths{vindex} = c(1:ds:end);
                else
                    paths{vindex} = c;
                end
                
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
end