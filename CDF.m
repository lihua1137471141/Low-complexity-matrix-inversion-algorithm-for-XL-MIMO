clear all;
clc;
%实验基本参数
%两个用户组
K1=16;
K2=16;
K=32;
%三个天线子阵列
M1=64;
M2=64;
M3=64;
S=3;
%M=64*S;
%Ms=64;
%三个VR的功率
p1=1;
P2=1;
P3=1;
p=1;
S=3;
sigma2 = 1;

numSetups=2;
snr = -5:5:30;
noiseVariancedBm = -200;
diagNorm='Norm1';
numRealizations=100;
Mrange=420:60:720;
a1=zeros(1,numRealizations);
a2=zeros(1,numRealizations);
a3=zeros(1,numRealizations);
a4=zeros(1,numRealizations);
time_elapsed = 0;
tic

   M=720;
        Ms=M/S;
        %Generate mobile communication setup
        [channelGaindB,R] = functionExampleSetup(M,S,K,diagNorm);
        %Compute channel gain over noise

        channelGainOverNoise = zeros(Ms,K,S);
        channelGainOverNoise(channelGaindB ~= 0) = channelGaindB(channelGaindB ~= 0) - noiseVariancedBm;

        %Go through each subarray

        H = functionChannelRealizations(Ms,K,channelGainOverNoise(:,:,1),R(:,:,:,1),numRealizations);
        H1=H(:,:,1:16);
        H = functionChannelRealizations(Ms,K,channelGainOverNoise(:,:,2),R(:,:,:,2),numRealizations);
        H2=H;
        H = functionChannelRealizations(Ms,K,channelGainOverNoise(:,:,3),R(:,:,:,3),numRealizations);
        H3=H(:,:,17:32);
        for n=1:numRealizations
            Hn1 = reshape(H1(:,n,:),[Ms K/2])';
            A=Hn1*Hn1'+ M/p*eye(K/2);
            X=[];
            for k=1:1:K/2
                e = zeros(K/2,1); e(k) = 1;
                x0=zeros(K/2,1);
                
                [x,B1]=JOR(A,e,x0,5,'',0.3);
                a_1= max(abs(eig(B1)));
                if a_1 >=a1(n)
                    a1(n)=a_1;
                end
                [x,B2]=JOR(A,e,x0,5,'',0.6);
                a_2= max(abs(eig(B2)));
                if a_2 >=a2(n)
                    a2(n)=a_2;
                end
                [x,B3]=JOR(A,e,x0,5,'',0.9);
                a_3= max(abs(eig(B3)));
                if a_3 >=a3(n)
                    a3(n)=a_3;
                end
                [x,B4]=JOR(A,e,x0,5,'',0.8);
                a_4= max(abs(eig(B4)));
                if a_4 >=a4(n)
                    a4(n)=a_4;
                end
            end
        end

    
    



[PDF1,Rate1] = hist( a1,40);
[PDF2,Rate2] = hist( a2,40);
[PDF3,Rate3] = hist( a3,40);
[PDF4,Rate4] = hist( a4,40);


PDF1 = PDF1/numRealizations;
PDF2 = PDF2/numRealizations;
PDF3 = PDF3/numRealizations;
PDF4 = PDF4/numRealizations;
for i=1:40
    CDF1(i) = sum(PDF1([1:i]));
end
for i=1:40
    CDF2(i) = sum(PDF2([1:i]));
end
for i=1:40
    CDF3(i) = sum(PDF3([1:i]));
end
for i=1:40
    CDF4(i) = sum(PDF4([1:i]));
end
figure
plot(Rate1,CDF1,'r-','LineWidth',2);
hold on

plot(Rate2,CDF2,'k-.','LineWidth',2);
%hold on
plot(Rate3,CDF3,'m--','LineWidth',2);
%hold on
%plot(Rate4,CDF4,'m--','LineWidth',2);
%     hold off
xlabel('Special radius');
ylabel('CDF');
legend('JOR,w=0.3','JOR,w=0.6','JOR,w=0.9','Location','best')
grid on