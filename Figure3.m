
%Empty workspace and close figures
close all;
clear;


%Select length of coherence block
tau_c = 200;

%Define range of number of BS antennas
Mrange = 10:10:100;

%Define range of number of UEs
Krange = 1:4:80;

%Set number of cells considered in the M-MMSE scheme
L = 9;
numIter=10;

%% Consider M=100 and varying K
K = 80;
M = max(Mrange);

%Compute number of samples for uplink data
tau_u = (tau_c-K);


%Compute complexity of receive combining
receiverProcessing = tau_u.*K*M;

            
            

                
                
           
                
                %Compute rKA computational complexity
                rKA=receiverProcessing+M*numIter + M;
                
           
            
per=5; 
for k1 = 1:length(Krange)
k=Krange(k1);
%Add complexity of computing combining matrix
%complexity_MMMSE = receiverProcessing + L*K*(M^2+M)/2 + M^2*K + (M^3-M)/3; %M-MMSE
%complexity_SMMSE = receiverProcessing + 3*M^2*K/2 + M*K/2 + (M^3-M)/3; %S-MMSE
%complexity_RZF = receiverProcessing + 3*K.^2*M/2 + 3*M*K/2 + (K.^3-K)/3; %S-MMSE2
complexity_ZF(k1) = 4*k^3+k-1+k;
complexity_JC(k1) = 2*k^2+(8*k^2-8*k)*per;
complexity_GS(k1) = 4*k^3-3*k^2+k+(8*k^2-8*k)*per;
complexity_JOR(k1) = 2*k^2+k+1+(8*k^2)*per;
%complexity_CG(k1) = (3*k^2+6*k)*per*4;
complexity_PCG(k1) = (8*k^2+46*k-6)*per+4*k^2+2*k;
complexity_CG(k1) = (8*k^2+46*k-6)*per;
end
%complexity_MR = receiverProcessing; %MR



%Plot the simulation results for M=100 and varying K



figure;
plot(Krange,complexity_ZF,'r->','LineWidth',2);
hold on
plot(Krange,complexity_CG,'b-o','LineWidth',2);
plot(Krange,complexity_PCG,'b-','LineWidth',2);
%plot(Krange,complexity_JC,'b:','LineWidth',1);
hold on
plot(Krange,complexity_GS,'k-.','LineWidth',2);
plot(Krange,complexity_JOR,'k-*','LineWidth',2);
%hold on

%hold on


%set(gca,'FontSize',12);
%plot(K([1 5:5:40]),complexity_MMMSE([1 5:5:40]),'rd','LineWidth',1);
%plot(K([1 5:5:40]),complexity_MR([1 5:5:40]),'bs','LineWidth',1);

xlabel('Number of UEs (K)','Interpreter','latex');
ylabel('Computation complex (Flops)','Interpreter','latex');
%set(gca,'YScale','log');
set(gca,'xLim',[10,40]);
grid on
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
legend('RZF','CG','Jac-PCG','GS','JOR','latex');


