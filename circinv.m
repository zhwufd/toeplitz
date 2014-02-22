function x=circinv(c,b)
% c is the first column of circulant matrix.
% b is right hand side vector.
e=length(c)*ifft(c);%eigenvalue of circulant matrix.
x=real(fft(ifft(b)./e));
end
