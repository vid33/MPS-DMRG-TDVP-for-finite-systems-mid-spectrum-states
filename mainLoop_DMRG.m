
%Create a stop button to exit main loop nicely
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
   
    [ Ar ] = transform_to_Ar( A, Dvec );
    if MID_SPECTRUM == false
        Al = DMRG_sweep_l_to_r( Ar, H1S, H2S, Dvec );
        Ar = DMRG_sweep_r_to_l( Al, H1S, H2S, Dvec );
    elseif MID_SPECTRUM == true
        Al = DMRG_sweep_l_to_r( Ar, H1S, H2S, Dvec, productState );
        Ar = DMRG_sweep_r_to_l( Al, H1S, H2S, Dvec, productState );
        tmp = overlap_A_with_productState( Al, productState );
     %   tmp.*conj(tmp)
    end
    

    A = Ar;
    
    currentStep = currentStep+1;
    
end