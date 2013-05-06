clear all;close all

n1=2;n2=2;n3=2;n4=2;
Q1=laplacian(n1)+0.01*speye(n1); Q2=laplacian(n2)+0.01*speye(n2); Q3=laplacian(n3)+0.01*speye(n3);
Q4=speye(n4);
%Q1=randn(n1); Q2=randn(n2); Q3=randn(n3);

Qseig={sparse(Q1) sparse(Q2) sparse(Q3) sparse(Q4)};
Qs={(Q1) (Q2) (Q3) (Q4)};


Qkron=kron(Q1,kron(Q2,Q3));

v=1:n1*n2*n3;
v=v(:);

tic
b1=kronmult(Qs,v);
toc

tic
b2=kronmult_eigen(Qseig,v);
toc

tic
b3=Qkron*v;
toc

norm(b3-b2)
norm(b3-b1)