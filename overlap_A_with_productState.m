function [ overlap ] = overlap_A_with_productState( A, productState )

    [~, N] = size(A);
    [~, d] = size(A{1});

    currentContraction = Contract({ A{1}, productState{1} }, {[1, -1], 1});
    for kk=2:N
        if kk < N
            currentContraction = Contract({currentContraction, A{kk}, productState{kk}}, {1, [1, 2, -1], 2});
        elseif kk== N
            currentContraction = Contract({currentContraction, A{kk}, productState{kk}}, {1, [1, 2], 2});
        end
            
    end

    overlap = currentContraction;
    
end

