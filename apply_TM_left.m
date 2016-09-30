function [ rho_out ] = apply_TM_left(A, rho_in, k1, k2 )
    %Apply transfer matrix to density matrix at k1 from site k1 to site k2;
    %from right to left
    for kk=k1:-1:k2
        if kk>1
            rho_in = Contract({A{kk}, conj(A{kk}), rho_in}, {[-1, 3, 1], [-2, 3, 2], [1, 2]});
        elseif kk == 1
            rho_in = Contract({A{kk}, conj(A{kk}), rho_in}, {[3, 1], [3, 2], [1, 2]});
        end
    end

    rho_out = rho_in;

end

