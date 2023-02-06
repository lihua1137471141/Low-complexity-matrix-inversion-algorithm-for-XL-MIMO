function [x,iter,time,res,resvec]=CG1(A,b,x,max_it,tol)
D=diag(diag(A));
U=D-triu(A);
L=D-tril(A);
w=1;
M=(D-w*L)*inv(D)*(D-w*U)/(w*(2-w));%计算ssor预处理矩阵M
[x,fl0,rr0,it0,rv0] = pcg(A,b,'',max_it,M);
x_true=A^(-1)*b;
% %L = ichol(A);
% %[x,fl2,rr2,it2,rv2] = pcg(A,b,'',max_it,L,L');
% %% 预优共轭梯度算法 (By X.Feng, Qiao.H, Zhang.L)
% %% To solve the equation Ax=b;
% % A : input matrix;  
% % b:right vector;
% % x : initial vector;
% % M : Preprocessing matrix; 
% % max_it : the maximum number of iterations;
% % tol : Accuracy;
% % time : CPU time;
% % iter : number of iterations at termination;
% % res : the norm of residual vector at termination;
% % resvec : the residual vector at termaination; 
%  
% D = diag(diag(A));
% L = -tril(A,-1);
% U = -triu(A,1);
%    % M = diag(diag(A)); % where use Jacobi Preprocessing matrix ; 
% %w=0.9;
%  %M=(2-w)^(-1)*(D/w+L)*(D/w)^(-1)*(D/w+L)';
%  %M=D;
%  M=(D+L)*D^(-1)*(D+L');
% tic;
% r = b - A * x;
% z = M \ r;
% p = z;
% rho = z' * r;
% mr = norm(r); % The L^2 norm of the vector r;
% iter = 0;
% while (iter < max_it)
%     iter = iter + 1; 
%     u = A * p;
%     alpha = rho / ( p' * u )';
%     x = x + alpha * p;
%     r = r - alpha * u;
%     z = M \ r;
%     rho1 = z' * r;
%     beta = rho1 / rho;
%     p = z + beta * p;
%     res = norm (r) / mr;
%     resvec = res;
%     if (res < tol )
%         break;
%     end
%     rho = rho1;
% end
% time = toc;
