function [ X ] = calculate_X_Al(A, rhor, Vl, Fl, H1S, H2S )
 
    [~, N] = size(A);
    X = cell(1, N);
    
    %2-site left
    for kk=1:N-1    
        if kk==1
            tmp1 = Contract({A{kk}, conj(Vl{kk})}, {[-1, -2], [-3, -4]});
        else
            tmp1 = Contract({A{kk}, conj(Vl{kk})}, {[1, -1, -2], [1, -3, -4]});
        end
        tmp2 = Contract({tmp1, H2S{kk+1}}, {[1, -1, 2, -3], [1,  -2, 2, -4]});
            
        if kk==N-1
            tmpA = Contract({conj(rhor{kk}^(-1/2)), conj(A{kk+1})}, {[-1, 1], [1, -2]}); 
            X{kk} = Contract({tmp2, A{kk+1}, tmpA}, {[1, 2, -1, 3], [1, 2], [-2, 3]});
        else
            tmpA = Contract({conj(rhor{kk}^(-1/2)), conj(A{kk+1})}, {[-1, 1], [1, -2, -3]}); 
            tmpB = Contract({A{kk+1}, rhor{kk+1}, tmpA}, {[-1, -2, 1], [1, 2], [-3, -4, 2]});
            X{kk} = Contract({tmp2, tmpB}, {[1, 2, -1, 3], [1, 2, -2, 3] });
        end
    end
    
    
    %2-site right
    for kk=1:(N-1)
        if kk==1
            tmp1 = Contract({A{kk}, conj(A{kk})}, {[-1, -2], [-3, -4]});
        else 
            tmp1 = Contract({A{kk}, conj(A{kk})}, {[1, -1, -2], [1, -3, -4]});
        end
        tmp2 = Contract({tmp1, H2S{kk+1}}, {[1, -1, 2, -3], [1, -2, 2, -4]});

 %       tmpA = Contract({conj(rhol{kk}^(-1/2)), conj(Vl{kk+1})}, {[1, -1], [1 -2, -3]});

        if kk==N-1
            tmp3 = Contract({tmp2, conj(Vl{kk+1}), A{kk+1}}, {[1, 2, 3, 4], [3, 4, -1], [1, 2]});
        else
            tmpB = Contract({A{kk+1}, rhor{kk+1}^(1/2)}, {[-1, -2, 1], [1, -3]});
            tmp3 = Contract({tmp2, tmpB, conj(Vl{kk+1})}, {[1, 2, 3, 4], [1, 2, -2], [3, 4, -1]});    
        end

        if kk==N-1
            X{kk+1} = tmp3;
        else
            X{kk+1} = X{kk+1} + tmp3;
        end      
    end
    
    %1-site
    for kk=1:N
        if kk==1
            tmp1 = Contract({A{kk}, rhor{kk}^(1/2)}, {[-1, 1], [1, -2]});
        else
            tmp1 = Contract({A{kk}, rhor{kk}^(1/2)}, {[-1, -2, 1], [1, -3]});
        end
        if kk==1
            tmp2 = Contract({tmp1, H1S{kk}}, {[1, -2], [1, -1]});
        else
            tmp2 = Contract({tmp1, H1S{kk}}, {[-1, 1, -3], [1, -2]});
        end
        if kk==1
            X{kk} = X{kk} +  Contract({tmp2, conj(Vl{kk})}, {[1, -2],[1, -1]});
        else
            X{kk} = X{kk} + Contract({tmp2, conj(Vl{kk})}, {[1, 2, -2], [1, 2, -1]}); 
        end
    end
    
    %F term
    for kk=1:N-1
       % tmpA = Contract({conj(rhol{kk}^(-1/2)), conj(Vl{kk+1})}, {[1, -1], [1 -2, -3]});
        if kk==N-1
             tmp1 = Contract({Fl{kk}, A{kk+1}}, {[1, -1], [1, -2]}); 
             X{kk+1} = X{kk+1} + Contract({tmp1, conj(Vl{kk+1})}, {[1, 2], [1, 2, -1]});
        else
            tmp1 = Contract({Fl{kk}, A{kk+1}, rhor{kk+1}^(1/2)}, {[1, -1], [1, -2, 2], [2, -3]});
            X{kk+1} = X{kk+1} + Contract({tmp1, conj(Vl{kk+1})}, {[1, 2, -2], [1, 2, -1]});
        end
%kk
%Contract({conj(Vl{kk+1}), tmp2}, {[1, 2, -1], [1, 2, -2]})
    end
    

end

