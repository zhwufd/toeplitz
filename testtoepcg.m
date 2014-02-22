% test file for solving the s.p.d toeplitz matrix by
% using preconditioned conjugate gradient method.
%%
% test para 1
clear;clc;close all;
thres=1;
n=2000;
tol=1e-3;
type='c';
%%
% test para 2
clear;clc;close all;
thres=2;
n=10000;
tol=1e-5;
type='c';
%%
% solving  toeplitz matrix.
% construct toeplitz matrix and right hand side b.
tic;
[ac,ar]=spdtoep(thres,n);
b=(rand(n,1)-0.5)*2*thres;
t1=toc;
%
tic;
[x,k]=toepcg(ac,ar,b,tol,type);
maint=toc;
% compute residual.
A=toeplitz(ac,ar);
resnorm=norm(A*x-b,2);
resmax=max(abs(A*x-b));
%%