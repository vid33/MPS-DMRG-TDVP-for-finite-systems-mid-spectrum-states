function [maxOverlap, I] = maxMPSOverlap(A, Dvec, currentPosition, Veff)

    [~, N] = size(A);
    [~, d] = size(A{1});
    [~, Heff_dim] = size(Veff);
    
    overlap_vec = zeros(1, Heff_dim);
        
    Ainit = A;
    
    for kk=1:Heff_dim
        Atmp = Veff(:, kk);
        if currentPosition == 1
            Atmp = reshape(Atmp, d, Dvec(currentPosition)); 
        elseif currentPosition == N
            Atmp = reshape(Atmp, Dvec(currentPosition-1), d);        
        else %away from endpoints
            Atmp = reshape(Atmp, Dvec(currentPosition-1), d, Dvec(currentPosition));
        end
        A{currentPosition} = conj(Atmp);
        overlap_vec(kk) = overlap_A_B(A, Ainit);
    end
    
%    format long e
%    overlap_vec.*conj(overlap_vec)
%    format short
    
    [maxOverlap, I] = max(real(overlap_vec.*conj(overlap_vec)));
    
end

