function [arr_smp] = rand_sample(n_pop, n_smp, arr_smp)
%Function receives:
    % n_pop:    pop size
    % n_smp:    Number of individuals to sample
    % arr_smp:  Array that will contain the sampled elements, in arr_smp(1)
    % an element in the population that will not be sampled is set
    for i=2:n_smp+1
        smp = randi(n_pop-i+1);
        m=1;
        for j=1:n_pop
            for k=1:i
                if k~=i && arr_smp(k)==j
                    break;
                end
            end
            if k==i
                if m==smp
                    break;
                else
                    m = m+1;
                end
            end
        end
        arr_smp(i)=j;
    end
end