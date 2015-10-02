function sigf=fourier(sig)
	dim_in = sum(size(sig)>1);
	switch dim_in
		case 1
			sigf=fft(sig);
		case 2
			sigf=fft2(sig);
	end
end