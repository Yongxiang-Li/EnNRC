function [period, Q] = Mul_NRC(signal, period)
    Q = nan(size(period));
    for j = period(:)'
        Q(j) = nrc(signal, j);
    end
end

function [ Q, m] = nrc( noiseSignal, N )
%NRC Summary of this function goes here
%   Detailed explanation goes here
    Y =  cut_signal(N,noiseSignal);
    m = size(Y,2);%  k = floor(m/n);
    Yhat = mean(Y,2)';   l = N*m;
    C = Yhat*Yhat'/length(Yhat);
    V = sum(noiseSignal(1:l).^2)/l - C;
    Q = C - V/(m-1);
end

function [ shapeSignal ] = cut_signal( l, signal)
%CUT Summary of this function goes here
%   Detailed explanation goes here
    n = length(signal);
    if l >= n
        shapeSignal = signal;
    else
        shapeSignal = reshape(signal(1:l*floor(n/l)), l, floor(n/l));
    end
end