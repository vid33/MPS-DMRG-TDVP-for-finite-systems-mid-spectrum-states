function [ Fl_1S, Fl_2S ] = calculate_Fl_Al( A, H1S, H2S, N )
    %Calculate effective Hamiltonian, effective == summed up to site k

    if nargin < 5
        [~, N] = size(A);
    end
        
    Fl_1S = cell(1,N);
    Fl_2S = cell(1,N);
    
    F1S_cell = cell(1, N);
    F2S_cell = cell(1, N);
    
    F1S_cell{1} = Contract({A{1}, H1S{1}, conj(A{1})}, {[1, -1], [1, 2], [2, -2]});

    %F2S_cell{1} = 0;
    
    tmp2S = Contract({A{1}, H2S{2}, conj(A{1})}, {[1, -1], [1, -2, 2, -4], [2, -3]});
    F2S_cell{2} = Contract({tmp2S, A{2}, conj(A{2})}, {[1, 2, 3, 4], [1, 2, -1], [3, 4, -2]});
   
    for kk=2:N
        tmp1 = Contract({A{kk}, conj(A{kk})}, {[1, -1, -2], [1, -3, -4]});
        F1S_cell{kk} = Contract({tmp1, H1S{kk}}, {[1, -1, 2, -2], [1, 2]});
    end
    
    for kk=2:N-1
        tmp1 = Contract({A{kk}, conj(A{kk})}, {[1, -1, -2], [1, -3, -4]});
        tmp2 = Contract({tmp1, H2S{kk+1}}, {[1, -1, 2, -3], [1,  -2, 2, -4]});
        F2S_cell{kk+1} = Contract({tmp2, A{kk+1}, conj(A{kk+1})}, {[1, 2, 3, 4], [1, 2, -1], [3, 4, -2]});
    end

    for kk=1:N
        Fl_1S{kk} = F1S_cell{kk};
        for pp=1:kk
            if kk-pp > 0
                Fl_1S{kk} = Fl_1S{kk} + apply_TM_right(A, F1S_cell{kk-pp}, kk-pp+1, kk );
            end
        end
    end 
    
    [~, D2] = size(A{1}); Fl_2S{1} = zeros(D2);
    for kk=2:N
        Fl_2S{kk} = F2S_cell{kk};
        for pp=2:kk
            if kk-pp > 0
                Fl_2S{kk} = Fl_2S{kk} + apply_TM_right(A, F2S_cell{kk-pp+1}, kk-pp+2, kk );
            end
        end
    end 
    
end

