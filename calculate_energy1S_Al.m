function [ energy1S ] = calculate_energy1S_Al( A, rhor , H1S)

    [~, N] = size(A);

    %endpoints
    tmp1 = Contract({A{1}, H1S{1}, conj(A{1})}, {[1, -1], [1, 2], [2, -2]});
    energy1S = Contract({tmp1, rhor{1}}, {[1, 2], [1,2]});
    
    tmp1 = Contract({A{N}, H1S{N}, conj(A{N})}, { [-1, 1], [ 1, 2], [-2, 2] });
    energy1S = energy1S+ trace(tmp1);
    
    %everything else
    for kk=2:N-1
        tmp1 = Contract({A{kk}, conj(A{kk})}, {[1, -1, -2], [1, -3, -4]});
        tmp2 = Contract({tmp1, H1S{kk}}, {[1, -1, 2, -2], [1, 2]});
        
        energy1S = energy1S + Contract({tmp2, rhor{kk}}, {[1, 2], [1, 2]});
    end

end

