    function [aic] = calc_AIC(RSS, N, K)
        if (N/K >= 40)
            aic = N * log(RSS / N) + 2* K ;
        else 
            aic = N * log(RSS / N) + 2*K + 2*K * (K+1) / (N-K-1);
        end
    end