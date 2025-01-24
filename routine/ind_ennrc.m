function [ X, limit ] = ind_ennrc( signal, C, maxPeriod, pr )
%IND_ENNRC calculate the statistics and control limit
%   suppose they are iid
Qyt = zeros(maxPeriod, 1);
for p = 1 : maxPeriod
    Qyt(p) = en_nrc(signal, p);
end

Z = Qyt ./ sqrt(diag(C));
X = sum(Z.^2)';
limit = chi2inv(1 - pr, maxPeriod);

end

