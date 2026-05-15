function [muV, stdV] = upd_PV(muV, stdV, winner, looser, Np)
    [D,~]    = size(winner);
    for i=1:D
        win = winner(i);
        los = looser(i);
        mu  = muV(i);

        if (mu > win) && (mu - win > 2-mu + win)
            win = win + 2;
        elseif (mu < win) && (win - mu > 2-win + mu)
            win = win - 2;
        end


        if (mu > los) && (mu - los > 2-mu + los)
                los = los + 2;
        elseif (mu < los) && (los - mu > 1-los + mu+1)
                los = los - 2;
        end

        muvi = mu + (1/Np)*(win-los);
        muV(i) = muvi;

        if(muV(i) < -1)
            muV(i) = muV(i)+2;
        elseif (muV(i) > 1)
            muV(i) = muV(i)-2;
        end

        temp = stdV(i)^2 + mu^2 - muvi^2 + (1/Np)*(win^2 - los^2);
        temp    = sqrt(abs(temp));
        stdV(i) = min([10, temp]);
    end
end