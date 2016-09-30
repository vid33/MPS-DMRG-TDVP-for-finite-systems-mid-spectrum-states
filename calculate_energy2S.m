function [ energy2S ] = calculate_energy2S( A, rhol, rhor , H2S)

    [~, N] = size(A);

    %endpoints
    tmp1 = Contract({A{1}, H2S{2}, conj(A{1})}, {[1, -1], [1, -2, 2, -4], [2, -3]});
    tmp2 = Contract({tmp1, A{2}, conj(A{2})}, {[1, 2, 3, 4], [1, 2, -1], [3, 4, -2]});
    
    energy2S = Contract({tmp2, rhor{2}}, {[1, 2], [1,2]});
    
    tmp1 = Contract({A{N}, H2S{N}, conj(A{N})}, { [-2, 1], [ -1, 1, -3, 2], [-4, 2] });
    tmp2 = Contract({tmp1, A{N-1}, conj(A{N-1})}, {[1 2, 3, 4], [-1, 1, 2], [-2, 3, 4]});
    
    energy2S = energy2S+ Contract({rhol{N-2}, tmp2}, {[1, 2], [1, 2]});
    
    %everything else
    for kk=2:N-2
        tmp1 = Contract({rhol{kk-1}, A{kk}, conj(A{kk})}, {[1, 2], [1, -1, -2], [2, -3, -4]});
        tmp2 = Contract({tmp1, H2S{kk+1}}, {[1, -1, 2, -3], [1,  -2, 2, -4]});
        tmp3 = Contract({tmp2, A{kk+1}, conj(A{kk+1})}, {[1, 2, 3, 4], [1, 2, -1], [3, 4, -2]});
        
        energy2S = energy2S + Contract({tmp3, rhor{kk+1}}, {[1, 2], [1, 2]});
    end

end

