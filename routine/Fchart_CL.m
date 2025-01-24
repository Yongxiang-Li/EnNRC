function [control_line] = Fchart_CL(m,pr)
% calculate the control limit for Fisher's test
% m is the No. of ordinates, pr is the probability of a larger value
       warning off
       if nargin == 1, pr = 0.0027; end
       G = fliplr(linspace(1/m, 50/m, 1000));
       for g = G
           k = floor(1/g);
           temp = zeros(k,1);
           for j = 1:k
               temp(j) = (-1)^(j-1)*nchoosek(m,j)*(1-j*g)^(m-1);
           end
           temp(isnan(temp)|isinf(temp)) = 0;
           if sum(temp)>pr
               break
           end
       end % figure; plot(fliplr(linspace(1/m, 100/m, 100))',result) 
       control_line = m*g;
end
