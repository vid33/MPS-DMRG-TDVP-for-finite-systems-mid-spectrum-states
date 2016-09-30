function [ deltaA ] = calculate_dAl(X, Vl, rhor)

    [~, N] = size(X);
    deltaA = cell(1, N);

    for kk=1:N
        if kk==1
            tmpA = Contract({X{kk}, (rhor{kk}^(-1/2))}, {[-1, 1], [1, -2]});
            deltaA{kk} = Contract({Vl{kk}, tmpA }, {[ -1, 1], [1, -2]});
        elseif kk>1
            tmpB2 = Contract({X{kk}, (rhor{kk}^(-1/2))}, {[-1, 1], [1, -2]});
            deltaA{kk} = Contract({Vl{kk}, tmpB2 }, {[ -1, -2,  1], [1, -3]});
        end
            
    end
        
end

