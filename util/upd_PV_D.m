function [muV, stdV] = upd_PV_D(muV, stdV, winner, looser, Np)
    D    = size(winner, 1);
    beta = 1/Np;
    for i=1:D
        win = winner(i);
        los = looser(i);
        mu  = muV(i);

        muvi = mu + beta*(win - los);
        muV(i) = muvi;

        if(muV(i) < -1)
            muV(i) = -1;
%             muV(i) = muV(i)+2;
        elseif (muV(i) > 1)
            muV(i) = 1;
%             muV(i) = muV(i)-2;
        end

        temp = (stdV(i)^2) + mu^2 - muvi^2 + beta*(win^2 - los^2);
        temp    = sqrt(abs(temp));
        stdV(i) = min([10, temp]);
    end
end