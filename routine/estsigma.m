function sigma = estsigma(signal, period)
l = length(signal);
Ve = zeros(size(period));
for N = period'
    signal(l+N-1) = 0;
    %     xsignal(l+N-1) = 0;
    Y = reshape(signal(1:ceil(l/N)*N), N, ceil(l/N));
    N2 = numel(Y)-l;
    N1 = N-N2;
    if N2>0
        Y1 = Y(1:N1,:);
        Y2 = Y(N1+1:end,1:end-1);
        [~, Ve1] = nrc(Y1);    [~, Ve2] = nrc(Y2);
        m = size(Y2,2);
        Ve(N) = (Ve1*(m+1)/m + Ve2*m/(m-1))/l;
    else
        [~, Ve(N)] = nrc(Y);
        m = size(Y,2);
        Ve(N) = Ve(N)*m/(m-1)/l;
    end
end
sigma = sqrt(min(Ve));
end