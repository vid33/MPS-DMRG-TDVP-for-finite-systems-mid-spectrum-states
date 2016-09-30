function [ Fr_1S, Fr_2S ] = calculate_Fr_Ar( A, H1S, H2S, N )
    %Calculate effective Hamiltonian, effective == summed up to site k

    if nargin < 5
        [~, N] = size(A);
    end
        
    Fr_1S = cell(1,N);
    Fr_2S = cell(1,N);
    
    F1S_cell = cell(1, N);
    F2S_cell = cell(1, N);
    
    F1S_cell{N} = Contract({A{N}, H1S{N}, conj(A{N})}, {[-1, 1], [1, 2], [-2, 2]});

    %F2S_cell{1} = 0;
    
    tmp2S = Contract({A{N}, H2S{N}, conj(A{N})}, {[-1, 1], [-2, 1, -4, 2], [-3, 2]});
    F2S_cell{N-1} = Contract({A{N-1}, conj(A{N-1}), tmp2S}, {[-1, 2, 1], [-2, 4, 3], [1, 2, 3, 4]});
   
    for kk=(N-1):-1:1
        if kk>1
            tmp1 = Contract({ A{kk}, conj(A{kk}) }, {[-2, -1, 1], [-4, -3, 1]});
            F1S_cell{kk} = Contract({H1S{kk}, tmp1}, {[1, 2], [1, -1, 2, -2] });
        elseif kk == 1
            tmp1 = Contract({ A{kk}, conj(A{kk})}, {[-1, 1], [ -2, 1]});
            F1S_cell{kk} = Contract({H1S{kk}, tmp1}, {[1, 2], [1, 2] });
        end
    end
    
    for kk=(N-2):-1:1
        tmp1 = Contract({A{kk+1}, conj(A{kk+1})}, {[-2, -1, 1], [-4, -3, 1]});
        tmp2 = Contract({H2S{kk+1}, tmp1}, {[-2, 1, -4, 2], [1,  -1, 2, -3]});
        if kk>1
            F2S_cell{kk} = Contract({A{kk}, conj(A{kk}), tmp2}, {[-1, 1, 2], [-2, 4, 3], [2, 1, 3, 4] });
        elseif kk == 1
            F2S_cell{kk} = Contract({A{kk}, conj(A{kk}), tmp2}, {[1, 2], [4, 3], [2, 1, 3, 4] });
        end
    end

    for kk=N:-1:1
        Fr_1S{kk} = F1S_cell{kk};
        for pp=N:-1:kk
            if pp-kk > 0
                Fr_1S{kk} = Fr_1S{kk} + apply_TM_left(A, F1S_cell{pp}, pp-1, kk);
            end
        end
    end 
    
    [D2, ~] = size(A{N}); Fr_2S{N} = zeros(D2);
    for kk=(N-1):-1:1
        Fr_2S{kk} = F2S_cell{kk};
        for pp=(N-1):-1:kk
            if pp-kk > 0
               Fr_2S{kk} = Fr_2S{kk} + apply_TM_left(A, F2S_cell{pp}, pp-1, kk );
            end
        end
    end 
    
end

