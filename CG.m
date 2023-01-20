function [X] = CG(A, b, TOL, ITE, initX)
% 调用格式: [rankA, rankB, N, X, ite, tol] = CNumbericCGMIteration(A, b, TOL, ITE, initX)
%           [rankA, rankB, N, X, ite, tol] = CNumbericCGMIteration(A, b, TOL, ITE)
%
% 作者: 王瑞
% 时间: 2015.10.28 16:54 - 19:21
% 版本: Version 1.0
%
% 任务: 共轭梯度法(Conjugate Gradient Method, CGM)迭代法求解线性方程组的解 Ax = b
%       适用于系数矩阵为对称阵的线性方程组(函数内包含矩阵对称化 b = A'*b, A = A'*A)
%       构建 x(k+1) = x(k) + alpha(k)*p(k)
%       等同 MATLAB 内置函数 cgs
%
% 输入: A = 系数矩阵, 方阵
%       b = 常系数向量, 行向量
%       ITE = 迭代次数上限
%       TOL = 解的精度(范数)
%       initX = 初始解
%
% 输出: rankA = 系数矩阵 A 的秩
%       rankB = 增广矩阵 B 的秩, 其中 B = [A|b]
%       N = 方程组未知量个数
%       X = 方程组的解向量
%       ite = 求解的迭代次数
%       tol = 实际误差
% 
if nargin == 5
    X = initX;
elseif nargin == 4
    X = zeros(size(b));
end
%% 解的判定
B = [A, b];
rankA = rank(A);
rankB = rank(B);
N = length(b);
ite = 0;
tol = inf;
if rankA ~= rankB
    disp('None: 矩阵 A 的秩不等于 [A b] 的秩， 无解！');
    return ;
elseif rankA ~= N
    disp(['Waring: 矩阵的秩小于' num2str(N) '[' num2str(rankA) ']''，存在无穷多解。']);
end
if rankA == 1
    disp('别闹，用手算的。');
    return;
end
if rankA == N
   % disp(['Good: 矩阵的秩为' num2str(N) '，存在唯一解。']);
end
% 系数矩阵对称，放大 A' 倍
if rank(A-A') ~= 0
    b = A'*b;
    A = A'*A;
end
% GCM 迭代
r = b - A*X;
while ite < ITE
    err = r'*r;
    ite = ite + 1;
    if ite == 1
        p = r;
    else
        beta = err / errold;
        p = r + beta*p;
    end
    Ap = A*p;
    alpha = err / ((Ap)'*p);
    X = X + alpha*p;
    r = r - alpha*Ap;
    errold = err;
    tol = norm(r) / norm(b);
    if tol < TOL
        %disp('Excatly: 求得解。');
        break;
    end
end
if ite > ITE
    disp(['Message: 在 ' num2str(ITE) ' 次迭代过程中，无法求得解，'...
        '建议增大总的迭代次数 或者 分析算法的收敛性。']);
end
end