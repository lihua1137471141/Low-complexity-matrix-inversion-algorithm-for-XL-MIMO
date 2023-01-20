function P = V_CG(H,p,M,K)
A=H*H'+ M/p*eye(K);
X=[];
for k=1:1:K
 e = zeros(K,1); e(k) = 1;
 x0=zeros(K,1);
[x]=CG(A,e,'',5,x0);
X=[X,x];
end
pre = H'*X;

P = sqrt(p/trace(pre*pre'))*pre;