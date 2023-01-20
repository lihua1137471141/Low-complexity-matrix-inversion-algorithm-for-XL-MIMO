clc
clear
%% 报错检验

% %当矩阵A奇异时迭代产生错误（这里只以Jacobi迭代为例）
% A = [1 1 ; 2 2];
% b = [1 1]';
% x0 = [1 1]';
% jacobi(A,b,x0)
% 
% %当不能通过行交换使对角线上元素都变成非零时，产生错误。
% A = [1 1 1;
%      1 1 0;
%      0 1 0];
% b = [1 1 1]';
% x0 = [1 1 1]';
% jacobi(A,b,x0)
% 
% %当迭代不收敛时，该迭代方法失效
% A = [1 2 ; 3 4];
% b = [1 1]';
% x0 = [1 1]';
% jacobi(A,b,x0)

%% 造一个次序相容的A，且满足Jacobi迭代矩阵J的特征值都为实数，J谱半径小于1.
%在这样的条件下，我们看看SOR方法最佳松弛因子实验值与理论值是否相符。
% 实验中，可以控制相同步数，选精度最高因子，也可以选定精度，看那个因子对应迭代
%的步数最少。实验选到的最优因子和理论最优因子有一定差距，原因在于：1、所采用的精度
%标准，并非真正的精度，这里的精度和解的逼近只有一个粗糙关系。2、不同因子对应的迭代
%步数可能会相同，相同步数对应的精度也有可能对相同。
clc
clear
A = [-3 1 0 1;1 -3 0 0; 0 0 -3 1; 1 0 1 -3];%A次序相容
n = size(A,2);
b = ones(n,1);
x0 = ones(n,1);

D = diag(diag(A));
J = eye(n) - D^-1*A;
if isreal(eig(J)) && (max(abs(eig(J))) < 1)
sep = 0:0.001:2;
w0 = sep(2:end-1);
m = size(w0,2);
Ts = zeros(1,m);
step_sets = zeros(1,m);
h = 0;
for w = w0 %采用不同因子用迭代法尝试，在相同迭代步数下选取精度最高。
[~,out] = SOR(A,b,x0,w,1.0e-4,15);
h = h+1;
Ts(h) = out;
[~,out] = SOR(A,b,x0,w,1.0e-4);
step_sets(h) = out;
end
inds = find(Ts==min(Ts));
ind = round((max(inds) + min(inds))/2);
%当多个精度相同时，取中间一个。
wb = w0(ind);%w 均匀分割迭代产生的最优松弛因子
wb0 = 2/( 1+sqrt( 1-max(abs(eig(J)))^2 ) );% 理论最优松弛因子
wb_error = wb - wb0;%理论与实际差距
[~,ind1]=min(abs(w0 - wb0));%最接近理论最优因子的实验因子的位置
step_error = ind - ind1;%位置步数差距
p = 500;
plot(w0(ind1-p:ind1+p),step_sets(ind1-p:ind1+p),'-mo',...
                'MarkerSize',3)%绘图直观
            title('同精度所需迭代步数随松弛因子变化图')
end


%% 我们来比较Jacobi迭代、Gauss_Seidel迭代和最优因子SOR迭代的速度。由于三
%%个算法每一步迭代所需的运算量相差不大，而程序中使用的内置函数不一而同，所以我们不
%%能通过比较时间，可通过比较步数来看看他们之间差距。程序中，我们通过不同迭代步数下
%%达到的实时精度情况，来直观地看三种方法的收敛性能。
clc
clear
close all
%%造一个主对角占优的正定矩阵，使之三种方法都收敛。
n = 10;
A = rand(n,n);
A = A'+A + eye(n)*2*n;
x0 = ones(n,1);
true_x = rand(n,1)*10;%真解
b = A*true_x;

% 通过上面介绍的方法，求SOR方法的最优松弛因子。
sep = 0:0.001:2;
w0 = sep(2:end-1);
m = size(w0,2);
Ts = zeros(1,m);
step_sets = zeros(1,m);
h = 0;
for w = w0
[~,out] = SOR(A,b,x0,w,1.0e-4,15);
h = h+1;
Ts(h) = out;
[~,out] = SOR(A,b,x0,w,1.0e-4);
step_sets(h) = out;
end
inds = find(Ts==min(Ts));
ind = round((max(inds) + min(inds))/2);
%当多个精度相同时，取中间一个。
wb = w0(ind);%w 均匀分割迭代产生的最优松弛因子

%不同迭代方式的一个比较
iterNum_Down = 1; %设置迭代步上下限
iterNum_Up = 10;
eps = 0.0001;
Ts_j = [];
Ts_G = [];
Ts_S = [];
Xs_j = [];
Xs_G = [];
Xs_S = [];
Xs_C = [];
Xs_C1 = [];
Xs_J = [];

iterNums = iterNum_Down : iterNum_Up;
for iterNum = iterNums
[x_j,out_j] = jacobi(A,b,x0,'',iterNum);
Ts_j(end+1) = out_j;
Xs_j(:,end+1) = x_j;
[x_G,out_G] = Gauss_Seidel(A,b,x0,'',iterNum);
Ts_G(end+1) = out_G;
Xs_G(:,end+1) = x_G;
[x_S,out_S] = SOR(A,b,x0,wb,'',iterNum);
Ts_S(end+1) = out_S;
Xs_S(:,end+1) = x_S;
[x_C] = CG(A, b, '', iterNum, x0);
%Ts_C(end+1) = out_C;
Xs_C(:,end+1) = x_C;
[x_C1] = CG1(A,b,x0,iterNum,'');
%Ts_C(end+1) = out_C;
Xs_C1(:,end+1) = x_C1;
[x_J,~]=JOR(A,b,x0,iterNum,0.5);
Xs_J(:,end+1) = x_J;

end


for k = 1:iterNum_Up - iterNum_Down + 1
error_j(k)= norm((Xs_j(:,k) - true_x), 2);
error_G(k)= norm((Xs_G(:,k) - true_x), 2);
error_S(k)= norm((Xs_S(:,k) - true_x), 2);
error_C(k)= norm((Xs_C(:,k) - true_x), 2);
error_C1(k)= norm((Xs_C1(:,k) - true_x), 2);
error_J(k)= norm((Xs_J(:,k) - true_x), 2);
figure;
end
error_m = norm(A\b - true_x);

plot(iterNums,error_J,'m-.','LineWidth',2)
 hold on
plot(iterNums,error_G,'c->','LineWidth',2)
%plot(iterNums,error_S,'b-.^','LineWidth',2)
%plot(iterNums,error_m,'b--o','LineWidth',2)
plot(iterNums,error_C,'bo','LineWidth',1)
plot(iterNums,error_C1,'k-','LineWidth',2)
grid on
legend('JOR','GS','CG','PCG');
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
xlabel('iteration','Interpreter','latex')
ylabel('Least square error','Interpreter','latex')
set(gca,'xLim',[1,6]);
%由于收敛速度都是呈高数量级的且不同方法收敛速度的巨大差异（特别是Jacobi迭代在图上
%降低了分辨率），所以我们很难从图上直观看出我们这三种方法随着迭代的进行，细致
%地看出收敛性能的不同。但是，我们可以从数值上简单地看到，通过十几步的迭代，
%三种方法的计算精度就已经达到了matlab数据类型（显示上的）精度（e-15左右）的上限。





function [x,out] = jacobi(A,b,x0,eps, iterNum)
%输入参数eps表示精度，iterNum表示迭代步数
%当输入参数为五个时，eps失效，以迭代步数作为终止条件
%% 
n = length(b);
if det(A) == 0 %借用了内置函数，行列式可用循环求解
    error('No Unique Solution Or No Solution!');
end
%% 当主对角线上元素不全为零时，可以通过交换两行来实现
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
            %若无论如何都无法通过交换两行使得对角线上元素为0时，原则上可以通过交换
            %两列实现，此时求得的结果对应元素也需交换。这里若通过交换行不能使得主对
            %角线非0，产生错误信息。
        end
    end
end
%% 针对不同的参数输入，来选择精度控制模式或者步数控制模式。
D = diag(diag(A));
invD = diag(1./diag(A));
J = invD * (D - A);
if max(abs(eig(J))) >= 1
    error('Jacobi Algorithm Cannot Convergence！')
end
%这里利用了Jacobi迭代收敛的充要条件判断是否收敛，求谱半径内置函数。事实上，可以
%通过矩阵分解、反幂法等手段求特征值。此判断并非必要，只是给A的选取是否合适，作为
%一个参考。收敛与否可通过最后x的变化来判断。
invDb = invD * b;
f = @(x)(J*x + invDb);
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
        eps = 1.0e-6;%缺省精度
    end
    out = 0;
    while 1
        x_0 = x;
        x = f(x);
        out = out + 1 ;
        if norm( x - x_0 ,2 ) < eps
            %这里以最后向量之间的距离来作为迭代精度，事实上不大精准。
            break;
        end
    end
end
end




function [x,out] = Gauss_Seidel(A,b,x0,eps,iterNum)
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


function [x,out] = SOR(A,b,x0,w,eps,iterNum)
%% 同Gauss_Seidel迭代，只是迭代函数f有所改变。增加了松弛因子w作为参数。
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
% if strcmp(w ,'best') %最佳松弛因子，对矩阵的要求较高。不在程序中考虑。
%     w = 2/(   1+sqrt(   1-max(abs(eig(J)))^2   )   );
% end
w_ = 1 - w;
invDb = invD * b;
L = invD * (D - tril(A));
B_w = D*(eye(n) - w*L)/w;
H_w = eye(n) - B_w^-1*A;
if max(abs(eig(H_w))) >= 1
    error('SOR Algorithm Cannot Convergence！')
end

    function x = f(x)
        for l = 1:n
            temp_x = J(l,:)*x + invDb(l);
            x(l) = w_*x(l) + w*temp_x;
        end
    end
x = x0;
if nargin == 6
    for k = 1:iterNum-1
        x = f(x);
    end
    x_0 = x;
    x = f(x);
    out = norm(x-x_0,2);
else
    if nargin == 4
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