
function [x,iter]=JOR(A,b,x0,nmax,xx,omega)
%JOR JOR method
% [X,ITER]=JOR(A,B,X0,NMAX,TOL,OMEGA) attempts to solve the system
% A*X=B with the JOR method. TOL specifies the tolerance of the method.
% NMAX specifies the maximum number of iterations. X0 specifies the initial
% guess. OMEGA is the relaxation parameter. ITER is the iteration number at
% which X is computed.
[n,m]=size(A);
if n ~= m, error('Only square systems'); end
iter=0;
r = b-A*x0; r0=norm(r); err=norm(r); x=x0;
while iter < nmax
iter = iter + 1;
for i=1:n
s = 0;
for j = 1:i-1, s=s+A(i,j)*x(j); end
for j = i+1:n, s=s+A(i,j)*x(j); end
xnew(i,1)=omega*(b(i)-s)/A(i,i)+(1-omega)*x(i);
end
x=xnew; r=b-A*x; err=norm(r)/r0;
end
return