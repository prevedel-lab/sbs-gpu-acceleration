function [bic] = calc_BIC(RSS, N, K)
    bic = N * log(RSS / N) + K * log(N);
end