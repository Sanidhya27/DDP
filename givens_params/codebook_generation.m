function codebook_generation( t, bits)
	x = rand(10000,1);
	c = [];

	for k=1:t-1
		y = asin(x.^(1/(2*k)));
		[tables,codes] = lloyds(y,2**bits);
		c = [c,codes];
	endfor
	f = strcat(num2str(bits),'bits.it');
	itsave(f,c);
end