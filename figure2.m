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
numSetups=5;
snr = 15:3:30;
noiseVariancedBm = -50;
diagNorm='Norm1';
numRealizations=100;
Mrange=420:60:720;

for s = 1:numSetups
    
        M=768;
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
     %   H2 = functionChannelRealizations(Ms,K,channelGainOverNoise(:,1:32,2),R(:,:,1:32,2),numRealizations);
        %H3=H1;
       % H3 = functionChannelRealizations(Ms,K,channelGainOverNoise(:,17:32,1),R(:,:,:,2),numRealizations);
      %  H3 = functionChannelRealizations(Ms,K/2,channelGainOverNoise(:,17:32,3),R(:,:,17:32,3),numRealizations);
      for m = 1:1:length(snr)
        Pow = sigma2*10^(snr(m)/10);
        sumSINR1(m,s) = RZF(Ms,K/2,Pow,numRealizations,H1,H2,H3,sigma2);
        sumSINR2(m,s) = RZF_CG(Ms,K/2,Pow,numRealizations,H1,H2,H3,sigma2);
        sumSINR3(m,s) = RZF_PCG(Ms,K/2,Pow,numRealizations,H1,H2,H3,sigma2);
        sumSINR4(m,s) = RZF_GS(Ms,K/2,Pow,numRealizations,H1,H2,H3,sigma2);
        sumSINR5(m,s) = RZF_JOR(Ms,K/2,Pow,numRealizations,H1,H2,H3,sigma2);
        

      end
    
end

figure
plot(snr,mean(sumSINR1,2)/3,'r->','LineWidth',2);
hold on
plot(snr,mean(sumSINR2,2)/3,'b-o','LineWidth',2);
plot(snr,mean(sumSINR3,2)/3,'b-','LineWidth',2);
plot(snr,mean(sumSINR4,2)/3,'k-.','LineWidth',2);
plot(snr,mean(sumSINR5,2)/3,'k-*','LineWidth',2);
set(gca,'xLim',[15,30]);
xlabel('SNR (dB)','Interpreter','latex')
ylabel('Average sum SE of each subarray (bps/Hz/subarray)','Interpreter','latex')
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
legend('RZF','CG (Non-stationary method)','SSOR-PCG (Non-stationary method)','GS (Stationary method)','JOR (Stationary method)','latex');
grid on