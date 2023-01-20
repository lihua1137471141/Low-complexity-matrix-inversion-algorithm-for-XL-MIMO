function [x,iter,time,res,resvec]=CG1(A,b,x,max_it,tol)
%% 预优共轭梯度算法 (By X.Feng, Qiao.H, Zhang.L)
%% To solve the equation Ax=b;
% A : input matrix;  
% b:right vector;
% x : initial vector;
% M : Preprocessing matrix; 
% max_it : the maximum number of iterations;
% tol : Accuracy;
% time : CPU time;
% iter : number of iterations at termination;
% res : the norm of residual vector at termination;
% resvec : the residual vector at termaination; 
 

    M = diag(diag(A)); % where use Jacobi Preprocessing matrix ; 

 
tic;
r = b - A * x;
z = M \ r;
p = z;
rho = z' * r;
mr = norm(r); % The L^2 norm of the vector r;
iter = 0;
while (iter < max_it)
    iter = iter + 1; 
    u = A * p;
    alpha = rho / ( p' * u )';
    x = x + alpha * p;
    r = r - alpha * u;
    z = M \ r;
    rho1 = z' * r;
    beta = rho1 / rho;
    p = z + beta * p;
    res = norm (r) / mr;
    resvec = res;
    if (res < tol )
        break;
    end
    rho = rho1;
end
time = toc;
