% NANOSCAT, (c) 2015 carmine e. cella
%
% pedagogical implementation of scattering transform
% (dyadic wavelets)

%% parameters

M = 2;					% Scattering order
J = 19;					% Maximal scale corresponding to T=Q*2^(J/Q+1)

%% load audio

[y, sr] = audioread ('../../datasets/various_data/5_notes.wav');
y = y(1:2^floor(log2(length(y)))); % truncate to power of 2
y = y / norm(y);

%% make filters (spline)

N = length (y);
a = 2;
spline_order=1;

S = @(omega)((1+2*cos(omega/2).^2)./(48*sin(omega/2).^4));
S2pi = 2^4;

psif = {};
phif = {};

% For all possible resolutions (original size divided by powers of 2)
for j0 = 0:floor(log2(N))
    N0 = N/2^j0;
    
    % Have we gone further than the largest scale?
    if N0 <= N/2^J;
        continue;
    end
    
    epsilon = 0;
    
    omega = [0:N0-1]'/N0*2*pi;
    
    for j1 = j0:J-1
        % Define spline wavelet
        if j1 == j0
            omega1 = a^(j1+1-j0)*omega;
            psif{j0+1}{j1+1}{1} = sqrt(2)*exp(-i*epsilon*omega1/2)./omega1.^(spline_order+1).* ...
                sqrt(S(omega1/2+pi)./(S(omega1).*S(omega1/2)));
            psif{j0+1}{j1+1}{1}(1) = 0;
            k2pi = N0/a^(j1+1-j0);
            pts = 1+k2pi:k2pi:N0/2+1;
            psif{j0+1}{j1+1}{1}(pts) = sqrt(2)*exp(-i*epsilon*omega1(pts)/2)./omega1(pts).^(spline_order+1).* ...
                sqrt(S2pi./S(omega1(pts)/2));
        else
            omega1 = a^(j1+1-j0)*omega(1:N0/2);
            psif{j0+1}{j1+1}{1} = [sqrt(2)*exp(-i*epsilon*omega1/2)./omega1.^(spline_order+1).* ...
                sqrt(S(omega1/2+pi)./(S(omega1).*S(omega1/2))) ; zeros(N0/2,1)];
            psif{j0+1}{j1+1}{1}(1) = 0;
            k2pi = N0/a^(j1+1-j0);
            pts = 1+k2pi:k2pi:N0/2;
            psif{j0+1}{j1+1}{1}(pts) = sqrt(2)*exp(-i*epsilon*omega1(pts)/2)./omega1(pts).^(spline_order+1).* ...
                sqrt(S2pi./S(omega1(pts)/2));
        end
    end
    
    % Define lowpass filter
    omega = [0:N0/2 -N0/2+1:-1]'/N0*2*pi;
    omega1 = a^(J-j0)*omega;
    phif{j0+1} = 1./(omega1.^(spline_order+1).*sqrt(S(omega1)));
    phif{j0+1}(1) = 1;
end

%% plot filters
close all
figure
for m = 1:numel(psif)
    figure
    for i = m:numel(psif{m})
        fprintf ('psi %d %d, length = %d\n', m, i, length (psif{m}{i}{1}));
        v = psif{m}{i}{1};
        plot (v)
        hold on
    end
end
%%
figure
for i = 1:numel(phif)
    plot (phif{i})
    fprintf ('phi %d, length = %d\n', i, length(phif{i}))
    hold on
end

%% convolve and subsample

aa = 1;
aa_psi = aa;
delta=1/log2(a);

next_bands=@(j) (max(0,j+delta));
downsampling_fac =@(res,j) max(0,floor(j * log2(a)-res-aa));
downsampling_fac_psi = @(res,j) max(0,floor(j * log2(a)-res-aa_psi));

%%
SC{1}{1}.signal = y;
SC{1}{1}.meta.scale=-1;
SC{1}{1}.meta.orientation=0;
SC{1}{1}.meta.resolution=0;

for m=1:M+1
    raster1=1;
    if m > size(SC,2)		% No more coefficients to calculate
        continue;
    end

    for s=1:numel(SC{m})
        sig=SC{m}{s}.signal;
        infos=SC{m}{s}.meta;
        
        % Decompose/propagate - retrieve high frequencies
        if m<=M            
            number_of_j = length(psif{1});
            sigf=fft(sig);
            res=infos.resolution;
            prev_j=(infos.scale>=0)*mod(infos.scale,number_of_j) + -100*(infos.scale<0);
            raster=1;
            children={};
            for j=max(0,next_bands(prev_j)):numel(psif{res+1})-1
                number_of_orientation = numel(psif{res+1}{j+1});
                for th=1:number_of_orientation
                    ds = downsampling_fac_psi(infos.resolution,j);
                    ds2 = 2^ds;
                    out = abs (ifft (sigf .* psif{res+1}{j+1}{th}));
                   if ds2 > 1
                       out = out(1:ds2:end)*sqrt(ds2);
                   end
                    
                    children{raster}.signal = out;
                    children{raster}.meta.resolution = res+ds;
                    children{raster}.meta.scale = (infos.scale>=0)*infos.scale*number_of_j + j;
                    children{raster}.meta.orientation = infos.orientation*number_of_orientation + th-1;
                    
                    raster=raster+1;
                end
            end
            for ic=1:numel(children)
                SC{m+1}{raster1}=children{ic};
                raster1=raster1+1;
            end
        end

        sigf=fft(sig);
        number_of_j = length(psif{1});
        smoothed_sig.meta=infos;
        ds = downsampling_fac(infos.resolution,number_of_j)^2;
        smoothed_sig.signal = ifft(sigf .* phif{infos.resolution+1});
       if ds > 1
           smoothed_sig.signal = smoothed_sig.signal(1:ds:end)*sqrt(ds);
       end
     
        smoothed_sig.meta.orignorm = norm(sig(:));
        SJ{m}{s} = smoothed_sig;
    
    end
end

%% format output -------------

yc = fft (y);
U1 = zeros(numel(psif{1}), length (y));
S1 = zeros(numel(psif{1}), length (y));

for i = 1:numel(psif{1})
    U1(i, :) = abs(ifft(yc .* psif{1}{i}{1}));
end

figure
imagesc (U1);

O1 = SC{2};
maxlen = length (O1{1}.signal);
nlambdas = numel(O1);

Scoeff = zeros(nlambdas, maxlen);

for i = 1:nlambdas
    itp = interpft(O1{i}.signal, maxlen);
    %sig = O1{i}.signal;
    %Scoeff(i, 1:length(sig)) = sig;
    Scoeff(i, :) = itp;
end

figure 
imagesc (Scoeff)

%%

for m = 1:M+1
    for i = 1:numel(SJ{m})
        fprintf ('order %d, lambda %d, coeff num = %d\n', m, i, length (SC{m}{i}.signal));
    end
end
