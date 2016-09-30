function [ A ] = transform_to_Al( A, Dvec )

    [~, N] = size(A);
    [~, d] = size(A{1});

    for kk=1:N
        if kk==1
            M = reshape(A{kk}, d, Dvec(kk));
        elseif kk == N
            M = reshape(A{kk}, d*Dvec(kk-1), 1);
        else
            M = reshape(A{kk}, d*Dvec(kk-1), Dvec(kk));
        end
            
        [Q, R] = qr(M, 0);
        
        if kk == 1
            A{kk} = reshape(Q, d, Dvec(kk));
        elseif kk == N
            A{kk} = reshape(Q, Dvec(kk-1), d); %NB rightmost R dropped so state normalised.
        else
            A{kk} = reshape(Q, Dvec(kk-1), d, Dvec(kk));
        end

        if kk < N
            A{kk+1} = Contract({R, A{kk+1}}, {[-1, 1], [1, -2, -3]});
        end
    end

end

