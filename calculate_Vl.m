function [ Vl ] = calculate_Vl(A, rhol, Dvec)
    % left gauge condition, calculaiton of relevant null-spaces

    [~, N] = size(A);
    [~, d, ~] = size(A{2});
    
    Vl = cell(1, N); %null space of contraction up to site kk...

    Vl{1} = zeros(2,2);

    for kk=2:N
        [rholD, ~] = size(rhol{kk-1});
    
        if kk< N 
            Dvec_next = Dvec(kk);
        elseif kk==N
            Dvec_next = 1; 
        end
  
        if d*rholD <= Dvec_next
            Vl{kk} = zeros(rholD, d, Dvec(kk));
        elseif d*rholD > Dvec_next
            sqrt_rhol_conjA = Contract({rhol{kk-1}^(1/2), conj(A{kk})}, {[-1, 1], [1, -2, -3]});
            sqrt_rhol_conjA = transpose(reshape(sqrt_rhol_conjA, d*rholD, Dvec_next ));
            Vl{kk} = null(sqrt_rhol_conjA);
            Vl{kk} = reshape(Vl{kk}, Dvec(kk-1), d, Dvec(kk-1)*d - Dvec_next);
          
            %test
            %sqrt_rhol_conjA = Contract({rhol{kk-1}^(1/2), conj(A{kk})}, {[-1, 1], [1, -2, -3]});
            %Contract({sqrt_rhol_conjA, Vl{kk}}, {[1, 2, -2], [1, 2, -1]})
        end
    end
 

end

