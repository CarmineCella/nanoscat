function smoothed_sig = smooth(sig,infos,phif, psif, downsampling_fac)
	% Smooth a signal and downsample
	
	sigf=fourier(sig);
	number_of_j = length(psif{1});
	smoothed_sig.meta=infos;
	ds = downsampling_fac(infos.resolution,number_of_j);
	smoothed_sig.signal = sub_conv(sig,sigf,phif{infos.resolution+1},2^ds);
	smoothed_sig.meta.orignorm = norm(sig(:));
end
