function x = sampleSolution(mu, sigma)
%sampleSolution This function is used to sample a decision variable on compact
%algorithms that use real representation
%   mu, is the median of the decision variable distribution.
%   sigma, is the standard deviation of the decision variable distribution.
%   Returns x, that is the sampled variable
%   This function is based on the code available in:
%   https://github.com/wkoder/mocde/tree/master
    if(sigma < 1e-8)
        x=mu;
    else
%         sqrt2Sigma  = sqrt(2) * sigma;
%         This is to use the integrated functions in matlab, some package
%         must be downloaded
%         erfMuNeg    = erf((mu-1)/sqrt2Sigma);
%         erfMuPlus   = erf((mu+1)/sqrt2Sigma);
% 
%         smp = rand();
%         C   = -erf((mu + 1) / sqrt2Sigma) / (erfMuNeg - erfMuPlus);
%         x   = mu - sqrt2Sigma * erfinv((smp - C) * (erfMuNeg - erfMuPlus));

        sqrt2Sigma  = 1.41421356237310 * sigma;
        invsqrt2Sigma = 1/sqrt2Sigma;
        erfMuPlus   = my_erf((mu+1)*invsqrt2Sigma);
        erfMuNeg    = my_erf((mu-1)*invsqrt2Sigma);
        smp = rand;
        if smp == 0
            smp = rand;
        elseif smp == 1
            smp = rand;
        end
        x   = mu - sqrt2Sigma * inv_erf(smp*(erfMuNeg - erfMuPlus) + erfMuPlus);
    end
end

