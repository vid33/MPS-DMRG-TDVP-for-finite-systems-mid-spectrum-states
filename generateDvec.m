function [ Dvec ] = generateDvec( d, Dmax, N )

    Dvec = zeros(1, N-1);
    
    for kk=1:N-1
        if kk < 15 % prevent matlab from calculating huge powers
            if d^(kk) < Dmax
                Dvec(kk) = d^(kk);
            else 
                Dvec(kk) = Dmax;
            end
        else 
            Dvec(kk) = Dmax;
        end
    end
    
    for kk=1:N-1
        if Dvec(N-kk) <= Dmax && kk < 15
            if Dvec(N-kk) > d^(kk)
                Dvec(N-kk) = d^(kk);
            end
        end
    end


end

