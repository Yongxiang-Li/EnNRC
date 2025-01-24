function [ control_limit ] = MaxNRC_CL( maxPeriod, pr )
%   Calculate the control limit for the Max NRC statistic
%   Suppose pr is small, which is ususally the case
    pr = 1 - pr;
    control_limit = norminv(pr);
    while normcdf(control_limit)^maxPeriod < pr
        control_limit = control_limit + 0.0001;
    end
end

