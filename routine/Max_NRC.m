function [ X, Z, Qyt ] = Max_NRC( signal, L, Periods )
%  Max_NRC calculate the maximum of the normalized feature
    signal = zscore(signal);
    Qyt = zeros(length(Periods), 1);
 
    for i = 1:length(Periods)
        Qyt(i) = en_nrc(signal, Periods(i));
    end

    Z = L \ Qyt;
    Z = zscore(flipud(Z));
    X = max(Z);
    % X = max(zscore(Z)); % the maximum of the feature, do zsocre

end

