function [x,H] = Gauss_Seidel(A,b,x0,eps,iterNum)
%% 同Jacobi迭代，只是迭代函数f有所改变。
n = length(b);
if det(A) == 0
    error('No Unique Solution Or No Solution!');
end
for i = 1 : n
    j = 1;
    while A(i,i) == 0
        if A(j,i) ~= 0 && A(i,j) ~= 0
            tempV = A(j,:);
            A(j,:) = A(i,:);
            A(i,:) = tempV;
            tempS = b(j);
            b(j) = b(i);
            b(i) = tempS;
        else
            j = j + 1;
        end
        if j > n
            error('No Exchange Exist!');
        end
    end
end
D = diag(diag(A));
invD = diag(1./diag(A));
J = invD * (D - A);
invDb = invD * b;
H = eye(n) - tril(A)^-1*A;
if max(abs(eig(H))) >= 1
    error('Gauss_Seidel Algorithm Cannot Convergence！')
end

    function x = f(x)
        for l = 1:n
            temp_x = J(l,:)*x + invDb(l);
            x(l) = temp_x;
        end
    end
x = x0;
if nargin == 5
    for k = 1:iterNum-1
        x = f(x);
    end
    x_0 = x;
    x = f(x);
    out = norm(x-x_0,2);
else
    if nargin == 3
        eps = 1.0e-6;
    end
    out = 0;
    while 1
        x_0 = x;
        x = f(x);
        out = out + 1 ;
        if norm( x - x_0 ,2 ) < eps
            break;
        end
    end
end
end