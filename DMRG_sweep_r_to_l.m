function [ Ar_out ] = DMRG_sweep_r_to_l( Al, H1S, H2S, Dvec, productState )

    [~, N] = size(Al);
    [d, ~] = size(Al{N});

    Ar = transform_to_Ar(Al, Dvec); %get rid of this later
    Al = transform_to_Al(Ar, Dvec); 
    
    A = Al;
        
   % rhor = calculate_rhor(Al, N);
    [ Fl_1S_Al, Fl_2S_Al ] = calculate_Fl_Al( A, H1S, H2S );
    Fl_Al = cell(1, N);
    Fl_Al{N-1} = Fl_1S_Al{N-1} + Fl_2S_Al{N-1};

    term1 =Contract({Fl_Al{N-1}, eye(d)}, {[-1, -3], [ -2, -4]});
    term1= reshape(term1, 4, 4);
    
    term2 = Contract({eye(Dvec(N-1)), H1S{N}}, {[-1, -3], [-2, -4]});
    term2 = reshape(term2, 4, 4);
    
    term3 = Contract({H2S{N}, A{N-1}, conj(A{N-1}), }, {[1, -2, 2, -4], [3, 1, -1], [3, 2, -3]});
    term3 = reshape(term3, 4, 4);

    Heff = term1 + term2 + term3;

    [Veff, Deff] = eig(Heff);
    if nargin < 5
        [~, I] = min(real(diag(Deff))); %'real' is important!
    else %MID_SPECTRUM == true
        [~, I] = maxMPSOverlap(A, Dvec, 1, Veff); 
    end
    Ac = Veff(:, I);
    Ac = reshape(Ac, Dvec(N-1), d);
    Ac = conj(Ac);
    
    M =   reshape(Ac, Dvec(N-1), d);
    [R, Q] = rq(M);
    A{N} = reshape(Q, Dvec(N-1), d);
    A{N-1} = Contract({A{N-1}, R}, {[-1, -2, 1], [1, -3]});
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %           [rhol] = calculate_rhol(A, N);
     %           [rhor] = calculate_rhor(A, N);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %           for zz=1:N+1
     %               test1(zz) = Contract({rhol{zz}, rhor{zz}}, {[1, 2], [1, 2]}); 
     %           end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %       energy1S = calculate_energy1S( A, rhol, rhor , H1S);
     %       energy2S = calculate_energy2S( A, rhol, rhor , H2S);
     %       energy1S+ energy2S
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %return;

            
    for kk=(N-1):-1:2
            [ Fr_1S_Ar, Fr_2S_Ar ] = calculate_Fr_Ar( A, H1S, H2S, kk );
            Fr_Ar_right = Fr_1S_Ar{kk+1} + Fr_2S_Ar{kk+1};
            
            term1 = Contract({eye(Dvec(kk-1)), eye(d), Fr_Ar_right }, { [-1, -4], [-2, -5], [-3, -6]});
            term1 = reshape(term1, Dvec(kk-1)*d*Dvec(kk), Dvec(kk-1)*d*Dvec(kk));
                
            Fl_Al{kk-1} = Fl_1S_Al{kk-1} + Fl_2S_Al{kk-1};
            term2 = Contract({Fl_Al{kk-1}, eye(d), eye(Dvec(kk)) }, { [-1, -4], [-2, -5], [-3, -6]});
            term2 = reshape(term2, Dvec(kk-1)*d*Dvec(kk), Dvec(kk-1)*d*Dvec(kk));
            
            if kk == N-1
                term3 = Contract({ eye(Dvec(kk-1)), H2S{kk+1}, A{kk+1}, conj(A{kk+1}) }, ...
                                            { [-1, -4 ], [-2, 1, -5, 2], [-3, 1], [-6, 2]}  );
            elseif  kk < N-1
                term3 = Contract({ eye(Dvec(kk-1)), H2S{kk+1}, A{kk+1}, conj(A{kk+1}) }, ...
                                            { [-1, -4 ], [-2, 1, -5, 2], [-3, 1, 3], [-6, 2, 3]}  );
            end
            term3 = reshape(term3, Dvec(kk-1)*d*Dvec(kk), Dvec(kk-1)*d*Dvec(kk));    
            
            if kk > 2
                term4 = Contract({ A{kk-1}, conj(A{kk-1}), H2S{kk}, eye(Dvec(kk)) }, ...
                                            { [1, 2, -1], [1, 3, -4], [2, -2, 3, -5], [-3, -6] });
            elseif kk == 2
                term4 = Contract({ A{kk-1}, conj(A{kk-1}), H2S{kk}, eye(Dvec(kk)) }, ...
                                            { [1, -1], [2, -4], [1, -2, 2, -5], [-3, -6] });       
            end
            term4 = reshape(term4, Dvec(kk-1)*d*Dvec(kk), Dvec(kk-1)*d*Dvec(kk));

            term5 = Contract({ eye(Dvec(kk-1)), H1S{kk}, eye(Dvec(kk)) }, { [-1, -4], [-2, -5], [-3, -6] } );
            term5 = reshape(term5, Dvec(kk-1)*d*Dvec(kk), Dvec(kk-1)*d*Dvec(kk));

            Heff = term1+ term2+ term3+ term4 + term5;
            
            [Veff, Deff] = eig(Heff);
            if nargin < 5
                [~, I] = min(real(diag(Deff))); %'real' is important!
            else %MID_SPECTRUM == true
                [~, I] = maxMPSOverlap(A, Dvec, kk, Veff); 
            end
            Ac = Veff(:, I);
            Ac = reshape(Ac, Dvec(kk-1), d, Dvec(kk));
            Ac = conj(Ac);            
            
            M =   reshape(Ac, Dvec(kk-1), d*Dvec(kk));
            [R, Q] = rq(M);
            A{kk} = reshape(Q, Dvec(kk-1), d, Dvec(kk));

            if kk > 2
                A{kk-1} = Contract({A{kk-1}, R}, {[-1, -2, 1], [1, -3]});
            elseif kk==2
                A{kk-1} = Contract({A{kk-1}, R}, {[-1, 1], [1, -2]});
            end
    end
        
    [ Fr_1S_Ar, Fr_2S_Ar ] = calculate_Fr_Ar( A, H1S, H2S, 1 );
    Fr_Ar_right = Fr_1S_Ar{2} + Fr_2S_Ar{2};

    term1 =Contract({eye(d), Fr_Ar_right}, {[-1, -3], [ -2, -4]});
    term1= reshape(term1, 4, 4); 

    term2 = Contract({H1S{1}, eye(Dvec(1)) }, {[-1, -3], [-2, -4]});
    term2 = reshape(term2, 4, 4);

    term3 = Contract({H2S{2}, A{2}, conj(A{2}) }, {[-1, 1, -3, 2], [-2, 1, 3], [-4, 2, 3]});
    term3 = reshape(term3, 4, 4);

    Heff = term1 + term2 + term3;

    [Veff, Deff] = eig(Heff);
    if nargin < 5
        [~, I] = min(real(diag(Deff))); %'real' is important!
    else %MID_SPECTRUM == true
        [~, I] = maxMPSOverlap(A, Dvec, N, Veff); 
    end
    Ac = Veff(:, I);
    Ac = reshape(Ac, Dvec(N-1), d);
    Ac = conj(Ac);

    M =   reshape(Ac, 1, d*Dvec(1));
    [R, Q] = rq(M);
    A{1} = reshape(Q, d, Dvec(1)); %drop R; it should be 1 anyway
                
    Ar_out = A;
    
end

