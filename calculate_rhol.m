%calculates left density matrix
function [ rhol ] = calculate_rhol( A, N )

    rhol = cell(1, N+1);  % rhol{N+1} := rhol{0}
    rhol{N+1} =1;

    rhol{1} = Contract({A{1}, conj(A{1})}, {[1, -1], [1, -2]});
    
    for kk=2:N-1
        rhol{kk} = Contract({rhol{kk-1}, A{kk}}, {[1, -1], [1, -2, -3]});
        rhol{kk} = Contract({rhol{kk}, conj(A{kk})}, {[1, 2, -1], [1, 2, -2]});    
    end
    
    rhol{N} = Contract({rhol{N-1}, A{N}}, {[1, -1], [1, -2]});
    rhol{N} = Contract({rhol{N}, conj(A{N})}, {[1, 2], [1, 2]});

end

