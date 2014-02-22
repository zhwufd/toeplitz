function [ac,ar]=spdtoep(thres,n)
% a generator of symmetric positive definite toeplitz matrix.
ac=zeros(n,1);
ac(2:n)=(rand(n-1,1)-0.5)*2*thres;
ac(1)=2*sum(abs(ac(2:n)))+50*rand(1,1);
ar=ac';
end
