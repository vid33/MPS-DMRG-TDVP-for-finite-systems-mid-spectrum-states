%calculates right density matrix
function [ rhor ] = calculate_rhor( A, N )

    rhor = cell(1, N+1);  % convention rhor{N+1} := rhor{0}
    rhor{N} =1;

    rhor{N-1} = Contract({A{N}, conj(A{N})}, {[-1, 1], [-2, 1]});  
    
    for kk=1:N-2
        rhor{N-kk-1} = Contract({A{N-kk}, rhor{N-kk}}, {[-1, -2, 1], [1, -3]});
        rhor{N-kk-1} = Contract({rhor{N-kk-1}, conj(A{N-kk})}, {[-1, 1, 2], [-2, 1, 2]});     
    end
    
    rhor{N+1} = Contract({A{1}, rhor{1}}, {[-1, 1], [1, -2]});
    rhor{N+1} = Contract({rhor{N+1}, conj(A{1})}, {[1, 2], [1, 2]});
    
end

