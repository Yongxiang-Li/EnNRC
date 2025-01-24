function [control_limit] = NRC_mcCL(len,Periods,pr,L)
%   Using the Monte Carlo method to estimate the control limit for NRC 

    pr = 1 - pr;
    N = 10000;
    cri_value = zeros(N,1);
    for i = 1:N
        signal = randn(len,1);
        signal = (signal - mean(signal))/std(signal);
        cri_value(i) = mul_ennrc( signal, L, Periods);
    end
    
    control_limit = quantile(cri_value,pr);

end