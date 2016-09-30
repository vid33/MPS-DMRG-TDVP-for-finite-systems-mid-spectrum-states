function [ A ] = generate_random_A_productState( Dvec, d )
    [~, N] = size(Dvec); N = N+1;
    A = cell(1,N);

    A{1} = zeros(d, Dvec(1));
    A{1}(1,1) = round(rand);
    if A{1}(1,1)==0
        A{1}(2,1) = 1;
    end
    A{N} = zeros(Dvec(N-1),d);
    A{N}(1,1) = round(rand);
    if A{N}(1,1) == 0
        A{N}(1,2) = 1;
    end
    for kk=2:(N-1)
        A{kk} = zeros(Dvec(kk-1),d,Dvec(kk));
        A{kk}(1,1,1) = round(rand);
        if A{kk}(1,1,1) == 0
            A{kk}(1,2,1) = 1;
        end
       % fprintf('AA\n');
       % A{kk}(1,1,1)
       % A{kk}(1,2,1)
    end
end

