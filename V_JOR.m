function P = V_JOR(H,p,M,K)
A=H*H'+ M/p*eye(K);
X=[];
for k=1:1:K
 e = zeros(K,1); e(k) = 1;
 x0=zeros(K,1);
[x]=JOR(A,e,x0,5,0.01,0.5);
X=[X,x];
end
pre = H'*X;

P = sqrt(p/trace(pre*pre'))*pre;