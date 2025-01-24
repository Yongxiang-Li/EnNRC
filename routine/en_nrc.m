function [ Q, C ] = en_nrc( noiseSignal, N )
%NRC Summary of this function goes here
%   Detailed explanation goes here

    l = length(noiseSignal);
    noiseSignal(l+N) = 0;

    m = floor(l/N);
    Y = reshape(noiseSignal(1:N*(m+1)), N, m+1); % cut signal
    N1 = mod(l,N); 
    Y1 = sum(Y(1:N1,:),2); 
    Y2 = sum(Y(N1+1:end,1:m),2); % cut into 2 parts
    C1 = Y1'*(Y1/l);
    C2 = Y2'*(Y2/l);
    C = C1/(m+1) + C2/m;
    V1 = sum(sum(Y(1:N1,:).^2)/l);
    V2 = sum(sum(Y(N1+1:end,1:m).^2)/l);
    Q1 = (C1-V1)/m;
    Q2 = (C2-V2)/(m-1);
    Q = (Q1 + Q2);
end