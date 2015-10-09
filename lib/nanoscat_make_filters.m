function [psi, phi, lp] = nanoscat_make_filters(N, J, shape)
if nargin < 3
    shape = 'gaussian';
end
nResolutions = 1 + floor(log2(N));
psi = cell(1, nResolutions);
phi = cell(1, nResolutions);

for res = 0:(nResolutions-1)
    N0 = N / 2^res;
    
    if N0 <= N/2^J;
        break;
    end
    
    for j = 0 : J -1
        sz = floor ((N0 /2^j) .* 0.8);
        if sz <= N/2^J
            break;
        end
        
        switch shape
            case 'hanning'
                v = zeros(N0, 1);
                v(1:(sz-1)) = (1 - cos(2*pi*(0:(sz-2))/(sz)));
            case 'gaussian'
                xi = 0.4 * 2^(-j);
                v = 2 * exp(- ((0:(N0-1))/N - xi).^2 * 10 * log(2) / xi^2 ).';
        end
        
        psi{res+1}{res+j+1} = v;
        
        if (res+j == J - 1)
            switch shape
                case 'hanning'
                    f = zeros(N0, 1);
                    half = floor (sz/2);
                    f(end-half+1:end) = v(1:half) * .5;
                    f(1:half)=v(half+1:half+half) * .5;
                    phi{res+1} = f;
                case 'gaussian'
                    bw = 0.4 * 2^(-1+J);
                    phi{res+1} = exp( - (0:(N0-1)).^2 * 10 * log(2) / bw^2).';
            end
        end
    end
end

lp = zeros (numel(psi{1}{1}), 1);
for i = 1 : numel(psi{1})
    lp = lp + 0.5 * (abs(psi{1}{i})).^2;
end
lp = lp + abs(phi{1}).^2;
lp(2:end) = (lp(2:end) + lp(end:-1:2)) * .5;
end
