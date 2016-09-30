function [ A ] = transform_to_Ar( A, Dvec )

    [~, N] = size(A);
    [~, d] = size(A{1});

    for kk=N:-1:1
        if kk==1
            M = reshape(A{kk}, 1, Dvec(kk)*d);
        elseif kk == N
            M = reshape(A{kk}, d, Dvec(kk-1));
        else
            M = reshape(A{kk}, Dvec(kk-1), d*Dvec(kk));
        end
            
        [R, Q] = rq(M);

        if kk == 1
            A{kk} = reshape(Q, d, Dvec(kk)); %NB leftmost R dropped so state normalised.
        elseif kk == N
            A{kk} = reshape(Q, Dvec(kk-1), d); 
        else
            A{kk} = reshape(Q, Dvec(kk-1), d, Dvec(kk));
        end

        if kk > 1
            if kk==2
                A{kk-1} = Contract({A{kk-1}, R}, {[-1, 1], [1, -2]});
            else
                A{kk-1} = Contract({A{kk-1}, R}, {[-1, -2, 1], [1, -3]});
            end
        end
    end

end

