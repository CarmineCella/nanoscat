function children = decompose(sig,infos,psif,next_bands,downsampling_fac)
	% Decompose a signal into its wavelet coefficients, downsample and compute the modulus
	
	number_of_j = length(psif{1});
	sigf=fourier(sig);
	res=infos.resolution;
	prev_j=(infos.scale>=0)*mod(infos.scale,number_of_j) + -100*(infos.scale<0);
	raster=1;
	children={};
	for j=max(0,next_bands(prev_j)):numel(psif{res+1})-1
		number_of_orientation = numel(psif{res+1}{j+1});
		for th=1:number_of_orientation
			ds = downsampling_fac(infos.resolution,j);
			out = abs(sub_conv(sig,sigf,psif{res+1}{j+1}{th},2^ds));

			children{raster}.signal = out;
			children{raster}.meta.resolution = res+ds;
			children{raster}.meta.scale = (infos.scale>=0)*infos.scale*number_of_j + j;
			children{raster}.meta.orientation = infos.orientation*number_of_orientation + th-1;

			raster=raster+1;
		end
	end
end