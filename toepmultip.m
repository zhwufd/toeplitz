function y=toepmultip(ac,ar,x)
n=length(x);
t=[ac;0;ar(end:-1:2)'];
xx=[x;zeros(n,1)];%pad zero to x.
e=2*n*ifft(t);%eigenvalue of circulant matrix.
yy=fft(e.*ifft(xx));
y=real(yy(1:n));
end