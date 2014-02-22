function [x,k]=toepcg(ac,ar,b,tol,type)
% use preconditioned conjugate gradient method to solve
%  symmetric positive definite toeplitz matrix.
%
% input:
% ac : the first column of toeplitz matrix.
% ar : the first row of toeplitz matrix.
% b : right hand side of toeplitz equation.
% tol : tolerance the accuracy of pcg method.
% use para type to choose which kind of precondition matrix.
% type can be the following option
%   's' : Strang preconditioned matrix.
%   'c' : T.Chan preconditioned matrix.
%   't' : Superoptimal preconditioned matrix.
%
% output:
% k : the number of iteration for preconditioned cg method.
% x : the solution of toeplitz equation.

% construct preconditioned matrix.
c=precondvec(ac,ar,type);% first column of precondition circulant matrix.
% pcg method.
x0=0;
r0=b;
p1=circinv(c,b);% p1=M\b
y0=p1;
k=0;
while norm(r0,2)>=tol
    z=toepmultip(ac,ar,p1);% z=A*p1
    v1=(y0'*r0)/(p1'*z);
    x0=x0+v1*p1;
    r1=r0-v1*z;
    y1=circinv(c,r1);% y1=M\r1
    mu=(y1'*r1)/(y0'*r0);
    p1=y1+mu*p1;
    y0=y1;
    r0=r1;
    k=k+1;
end
%
x=x0;