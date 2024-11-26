function [ M ] = get_cov_matrix( L, P )
% The covariance matrix of EnNRC
    len = length(P);
    M = nan(len, len);
    for i = 1:len
        p = P(i);
        op = ceil(L/p);    ep = floor(L/p);
        p1 = L-ep*p;    p2 = p-p1;
        for j = 1:len
            q = P(j);
            if ~isnan(M(i,j)), continue; end
            if mod(q, p)==0  % mod(q, p)==0
                M(i,j) = 2*(op*p1/(op-1)+ep*p2/(ep-1))/L^2;
            else
                oq = ceil(L/q);    eq = floor(L/q);
                q1 = L-eq*q;    q2 = q-q1;
                n = lcm(p,q); % n is the least common multiple of p and q
                if n <= L
                    on = ceil(L/n);    en = floor(L/n);    n1 = L-en*n;
                    Mnp = repmat([ones(p1,1)/(op-1); ones(p2,1)/(ep-1)], n/p, 1);
                    Mnq = repmat([ones(q1,1)/(oq-1); ones(q2,1)/(eq-1)], n/q, 1);
                    M(i,j) = 2*(en*(en-1)*(Mnp'*Mnq) + 2*en*(Mnp(1:n1)'*Mnq(1:n1)))/L^2;
                else
                    M(i,j) = 0;
                end
            end
            M(j,i) = M(i,j);
        end
    end
end

