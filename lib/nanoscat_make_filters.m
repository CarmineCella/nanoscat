function [psi, phi, lp] = nanoscat_make_filters(N, J)
nResolutions = 1 + floor(log2(N));
psi = cell(1, nResolutions);
phi = cell(1, nResolutions);

for res = 0:(nResolutions-1)
    N0 = N / 2^res;
    
    if N0 <= N/2^J;
        break;
    end
    
    for j = 0 : J -1
        v = zeros(1, N0);
        
        sz = floor ((N0 /2^j));
        if sz <= N/2^J;
            break;
        end
        
        s = (1 - cos(2*pi*(0:sz-1)'/(sz))); % hanning-zero-zero
        v(1:sz) = s;
        v(sz) = 0;
        
        psi{res+1}{res+j+1} = v';
        
        if (res+j == J - 1)
            f = zeros(1, N0);
            half = floor (sz/2);
            f(end-half+1:end) = v(1:half);
            f(1:half)=v(half+1:sz);
            phi{res+1} = f';
        end
    end
end

lp = zeros (numel(psi{1}{1}), 1);
for i = 1 : numel(psi{1})
    lp = lp + 0.5 * (abs(psi{1}{i})).^2;
end
lp = lp + abs(phi{1}).^2;
lp = (lp + lp(end:-1:1)) * .5;
end
