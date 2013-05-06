clear all;close all;

nx=64;
nshifts=9;
distance=1;
colouringtype=1;
verbosity=4;
cgTol=1e-6;
cgMaxIter=10000;

Q=modLaplacian([nx nx],'stat3');
Q=Q+(1e-2)*speye(size(Q));

trueMargVars=diag(inv(full(Q)));


Qtemp=Q^5;
probe=matlab_wrap_probing_vector(Qtemp,colouringtype,distance,verbosity);

d=eig(full(Q));
minEig=d(1);
maxEig=d(length(Q));

clear d;

[wsq,dzdt,intConst]=matlab_wrap_compute_shifts(minEig,maxEig,nshifts,verbosity);

ld=2*sum(log(diag(chol(Q))));


% ld1=matlab_wrap_eigen_det_app(Q,cgTol,cgMaxIter,distance,3,3);

% Q=triu(Q);

ld2=matlab_wrap_det_app_base(Q,probe,wsq,dzdt,intConst,cgTol,cgMaxIter,4);



Qtemp=Q^34;
probe=matlab_wrap_probing_vector(Qtemp,colouringtype,distance,3);


M=speye(size(Q));

margVars=matlab_wrap_margvars(Q,M,probe,cgTol,cgMaxIter);

relDiffDet=ld/ld2

plot(margVars./trueMargVars);


figure;plot(trueMargVars)
figure;plot(margVars)