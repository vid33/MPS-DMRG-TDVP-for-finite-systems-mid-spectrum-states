function [ A ] = generate_random_A( Dvec, d )
    [~, N] = size(Dvec); N = N+1;
    A = cell(1,N);

    A{1} = 0.5*(randn(d,Dvec(1))+1i*randn(d,Dvec(1)))/sqrt(d);
    A{N} = 0.5*(randn(Dvec(N-1),d)+1i*randn(Dvec(N-1),d))/sqrt(d);
    for kk=2:(N-1)
        Dav = (1/2)*(Dvec(kk) + Dvec(kk-1)); %this rescaling important, in order to have norm approx 1. Danger of blowup or decay exp in N.
        A{kk} = 0.5*(randn(Dvec(kk-1),d,Dvec(kk))+1i*randn(Dvec(kk-1),d,Dvec(kk)))/sqrt(Dav);
    end
end

