function [ Y ] = randGP(N, theta)
    %RANDQPGP generate QPGP signal
    r0 = gauss_cov(theta,(1:N)', 1);
    r0 = [r0; flipud(r0(2:end))];    r0(1) = r0(1) + 1e-10;
    S = fft(r0); % S = real(fft(r0));    S(S<0) = 0;
    Z = [randn(N,1); zeros(N-1,1)];
    Y = ifft(fft(Z).*sqrt(S));
    Y = Y(1:N);
    if nargout == 0,    figure; plot(Y);    end
%     C = toeplitz(r0);    D = dftmtx(length(r0));
%     Y = D'*((D*Z).*sqrt(S))/length(r0);    figure; plot(real(Y(1:N))); 
%     C = gauss_cov(theta,(1:N)', (1:N)');
%     L = chol(C+1e-10*eye(N), 'lower');
%     Y = L*Z(1:N);    figure; plot(Y);
end

function [ R ] = gauss_cov( theta, X, Y )
% p: is the period
% 
    if nargin == 3
        D = X - Y';
    else
        Y = X;
        D = X - Y';
    end
    D = theta * D;
    R = exp(-D.^2);
end







