function [ stat, I ] = Fchart( signal )

% Calculate the statistic of Fisher's test
% m is used to calculate No. of ordinates
    n = length(signal);
    m = ceil((n+1)/2);
%     A = signal'*exp(-1j*2*pi*((0:n-1)'*(1:m))/n);
%     I = (A.*conj(A))/n; % normalizition parameter is not important
    power = fft(signal); % this is even faster
    I = conj(power(2:m)).*power(2:m)/n;
    stat = max(I)/mean(I);
end
