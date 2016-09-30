function [ rho_out ] = apply_TM_right(A, rho_in, k1, k2 )
    %Apply transfer matrix to density matrix at k1-1 from site k1 to site k2

    for kk=k1:k2
        rho_in = Contract({rho_in, A{kk}, conj(A{kk})}, {[1, 2], [1, 3, -1], [2, 3, -2]});
    end

    rho_out = rho_in;

end

