function [ deltaA ] = calculate_dA(X, Vl, rhol, rhor)

    [~, N] = size(X);
    deltaA = cell(1, N);

    for kk=1:N
        if kk==1
            tmpA = Contract({X{kk}, (rhor{kk}^(-1/2))}, {[-1, 1], [1, -2]});
            deltaA{kk} = Contract({Vl{kk}, tmpA }, {[ -1, 1], [1, -2]});
        elseif kk>1
            tmpB1 = Contract({rhol{kk-1}^(-1/2), Vl{kk}}, {[1, -1], [1, -2, -3]});
            tmpB2 = Contract({X{kk}, (rhor{kk}^(-1/2))}, {[-1, 1], [1, -2]});
            deltaA{kk} = Contract({tmpB1, tmpB2 }, {[ -1, -2,  1], [1, -3]});
        end
            
    end
        
end

