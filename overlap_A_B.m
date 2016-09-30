function [ overlap ] = overlap_A_B( A, B )

    [~, N] = size(A);

    currentContraction = Contract({ A{1}, conj(B{1}) }, {[1, -1], [1, -2]});
    for kk=2:N-1
            currentContraction = Contract({currentContraction, A{kk}, conj(B{kk})}, {[1, 2], [1, 3, -1], [2, 3, -2]});            
    end

    currentContraction = Contract({currentContraction, A{N}, conj(B{N})}, {[1, 2], [1, 3], [2, 3]});
    
    overlap = currentContraction;
    
end

