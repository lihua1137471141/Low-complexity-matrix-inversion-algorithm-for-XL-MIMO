function sumSINR = RZF(M,K,p,numRealizations,H1,H2,H3,sigma2)
% Computes the sum of the signal-to-interference-plus-noise-ratio (SINR) of each user.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of BS antennas.
% @param  K               number of users.
% @param  p               total uplink transmit power per UE [mW].
% @param  numRealizations number of channel realizations (small-fading).
% @param  H               M x numRealizations x K matrix with channel responses.
% @return V               M x numRealizations x K RZF receive combining matrix.
%

%% Preamble

%Hold the K x K identity matrix
eyeK = eye(K);

%Prepare to save average sum SINR
sumSINR = zeros(K,numRealizations);

%Go through all channel realizations
for n = 1:numRealizations

    %Extract channel realizations from all users to the BS
    Hn1 = reshape(H1(:,n,:),[M K])';
    Hn2 = reshape(H2(:,n,:),[M 2*K])';
    Hn3 = reshape(H3(:,n,:),[M K])';
    %  Hn4 = reshape(H4(:,n,:),[M K]);
% Hn1 = sqrt(1/2)*(randn(K,M) + 1i*randn(K,M));
%    Hn2 = sqrt(1/2)*(randn(2*K,M) + 1i*randn(2*K,M));
%    Hn3 = sqrt(1/2)*(randn(K,M) + 1i*randn(K,M));


    %Compute RZF combining
    pre1 = Hn1'*inv(Hn1*Hn1'+ M/p*eye(K));
    V_RZF1 = sqrt(p/trace(pre1*pre1'))*pre1;
     pre2=Hn2'*inv(Hn2*Hn2'+M/p*eye(2*K));
    V_RZF2 =sqrt(p/trace(pre2*pre2'))*pre2;
     pre3=Hn3'*inv(Hn3*Hn3'+M/p*eye(K));
    V_RZF3 =sqrt(p/trace(pre3*pre3'))*pre3;
    % V_RZF4 = p*Hn4/(p*(Hn4'*Hn4)+eyeK);

    %Go through all users
    for k1 = 1:K
        signal1 = abs(Hn1(k1,:)*V_RZF1(:,k1)).^2;
        signal2 = abs(Hn2(k1,:)*V_RZF2(:,k1)).^2;
        %combiningNorm1 = norm(V_RZF1).^2+norm(V_RZF2).^2;
        interf1=0;
        interf2=0;
        for k2 = 1:K
        %Compute signal and interference + noise terms       
        interf1 = interf1 + abs(Hn1(k1,:)*V_RZF1(:,k2)).^2;
        end
        for k3 = 1:2*K
        interf2 = interf2 + abs(Hn2(k1,:)*V_RZF2(:,k3)).^2;
        end
        %Compute sum SINR
        sumSINR1(k1,n) =log2(1+(signal1+signal2 / (interf1+interf2 - signal1- signal2 + sigma2)));
        
    end



     for k1 = 1:K
        signal1 = abs(Hn3(k1,:)*V_RZF3(:,k1)).^2;
        signal2 = abs(Hn2(k1+8,:)*V_RZF2(:,k1+8)).^2;
        %combiningNorm1 = norm(V_RZF1).^2+norm(V_RZF2).^2;
        interf1=0;
        interf2=0;
        for k2 = 1:K
        %Compute signal and interference + noise terms       
        interf1 = interf1 + abs(Hn3(k1,:)*V_RZF3(:,k2)).^2;
        end
        for k3 = 1:2*K
        interf2 = interf2 + abs(Hn2(k1,:)*V_RZF2(:,k3)).^2;
        end
        %Compute sum SINR
        sumSINR2(k1,n) = log2(1+(signal1+signal2 / (interf1+interf2 - signal1- signal2 + sigma2)));
        
    end

end
%Treating NaNs
sumSINR1(isnan(sumSINR1)) = 0;
sumSINR2(isnan(sumSINR2)) = 0;

%Averaging and summing
sumSINR = sum(mean(sumSINR1,2))+sum(mean(sumSINR2,2));