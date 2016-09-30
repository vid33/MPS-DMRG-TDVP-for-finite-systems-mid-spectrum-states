
    testLGC = zeros(1,N);
    for kk=1:N
        if kk==1
            tmp = Contract({deltaA{kk}, conj(A{kk})}, {[1, -1], [1, -2]});
        else
            tmp = Contract({rhol{kk-1}, deltaA{kk}, conj(A{kk})}, {[1, 2], [1, 3, -1], [2, 3, -2]});
        end
        testLGC(kk) = Contract({tmp, rhor{kk}}, {[1,2], [1,2]});
    end

    norm_dA_density_test = zeros(1, N);

    for kk=1:N
        if kk==1
            tmp = Contract({deltaA{kk}, conj(deltaA{kk})}, {[1, -1], [1, -2]});
        else
            tmp = Contract({rhol{kk-1}, deltaA{kk}, conj(deltaA{kk})}, {[1, 2], [1, 3, -1], [2, 3, -2]});
        end
        norm_dA_density_test(kk) = Contract({tmp, rhor{kk}}, {[1,2], [1,2]});
    end

    clearvars tmp;


