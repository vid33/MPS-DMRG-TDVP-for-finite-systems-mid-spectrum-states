
%Create a stop button so loop below exits nicely
figh = figure;

global IS_ABORTED;
IS_ABORTED = false;

btn = uicontrol('style', 'pushb', 'string', 'Abort', ...
                'callback', @doAbort);
drawnow;

currentStep = 0;

while 1

    %left and right density matrices
    [ A ] = transform_to_Al( A, Dvec );

%    [rhol] = calculate_rhol(A, N);
    [rhor] = calculate_rhor(A, N);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk=1:N+1
       test1(kk) = trace(rhor{kk}); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    energy1S = calculate_energy1S_Al( A, rhor, H1S);
    energy2S = calculate_energy2S_Al( A, rhor, H2S);

    Vl= calculate_Vl_Al( A, Dvec );
    
    [ Fl_1S, Fl_2S ] = calculate_Fl_Al( A, H1S, H2S );

    for kk=1:N
        Fl{kk} = Fl_1S{kk} + Fl_2S{kk};
    end
    
    [ X ] = calculate_X_Al(A, rhor, Vl, Fl, H1S, H2S );
 
    for kk=1:N
        norm_dA_density(kk) = Contract({X{kk}, conj(X{kk})}, {[1,2], [1,2]});
    end

    [dA ] = calculate_dAl(X, Vl, rhor);

    norm_dA = sum(norm_dA_density);

    info;

    %exit nicely if button pressed
    drawnow;
    if IS_ABORTED == true
        close(figh);        
        break;
    end

    if rem(currentStep,250)==0
        save('data/tmp_save.mat');
    end 

    for kk=1:N
       A{kk} = A{kk} - dt*dA{kk};
    end
    
    currentStep = currentStep+1;
    
end