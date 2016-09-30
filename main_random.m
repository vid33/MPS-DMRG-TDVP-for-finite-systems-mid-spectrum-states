clear;

TDVP = false; %else DMRG
if TDVP == true
    MID_SPECTRUM = false;
end

%model set only one to true;
ISING = false;
HEISENBERG = false;
THERMALISATION_TEST = true;

N = 20; %no of sites
Dmax = 4; %max bond dimension
d = 2; %spin dimension

if TDVP == true
    dt = 0.01;
elseif TDVP == false %DMRG == true
    MID_SPECTRUM = true; %otherwise optimise for groundstates
    if MID_SPECTRUM == true
        productState = cell(1,N); %max overlap with this product state, rathern than min(Heff)
        for kk = 1:N
            rand_tmp = round(rand);
            if rand_tmp == 1  
                productState{kk} = [0; 1];
            else
                productState{kk} = [1; 0];
            end
        end
    end
end
    
sx=[0 1;1 0]; Sx = (1/2)*sx;
sy=[0 -1i; 1i 0]; Sy = (1/2)*sy;
sz=[1 0;0 -1]; Sz = (1/2)*sz;

H2S = cell(1, N); %two-site contribution to H
H1S = cell(1, N); %one-site contribution to H
Fl = cell(1,N);

norm_dA_density = zeros(1, N);

% H  = sum -Sz \otimes Sz    
% Exact energy -(1/4)*(N-1)*J
% Note D=2 saturates

if ISING == true
    J = ones(1, N); J(1) = 0;
    %h = 0.5*ones(1, N);
    h = 0.5- rand(1,N);
    for kk=1:N
        H2S{kk} = -J(kk)*Contract({Sz,Sz},{[-1,-3],[-2,-4]});
        H1S{kk} = -h(kk)*Sx;
    end
    fOut = sprintf('data/Ising/Dmax=%d/Dmax=%dN=%dJ=%dh=%d.mat', Dmax, Dmax, N, J(2), h(2));
elseif HEISENBERG == true
    %Exact energy for Jx=Jy=Jz=1, h=0  E= 14?ln(2)=?0.443147
    Jx = ones(1, N); Jx(1) = 0;
    Jy = ones(1, N); Jy(1) = 0;
    Jz = ones(1, N); Jz(1) = 0;
    %h = zeros(1, N);
    h = 8*2*(0.5- rand(1,N));
    for kk=1:N
        H2S{kk} = Jx(kk)*Contract({Sx,Sx},{[-1,-3],[-2,-4]}) ...
                    +Jy(kk)*Contract({Sy,Sy},{[-1,-3],[-2,-4]}) ...
                    +Jz(kk)*Contract({Sz,Sz},{[-1,-3],[-2,-4]});
        H1S{kk} = -h(kk)*Sz;
    end
    fOut = sprintf('data/Heisenberg/Dmax=%d/Dmax=%dN=%dJx=%dJy=%dJz=%dh=%d.mat', Dmax, Dmax, N, Jx(2), Jy(2), Jz(2), h(2));
elseif THERMALISATION_TEST == true  %% Philip's hamiltonian
    tt=-1*ones(1,N);  tt(1)=0;
    h=zeros(1,N);
    dh_sigma=5; %standard deviation of normal distribution
    for kk=1:N
       H2S{kk} = tt(kk)*Contract({(1/2)*(sz+eye(2)), sx}, {[-1,-3],[-2,-4]});
       H1S{kk} = ( h(kk)+normrnd(0, dh_sigma) )*(1/2)*(sz+eye(2));
    end
    
end

%fOut = sprintf('data/D=%d/d=%dD=%dh=%d.mat', Dmax, d, Dmax, h);
Dvec = generateDvec(d, Dmax, N); %specifies bond dimensions along chain
if MID_SPECTRUM == true
    [A] = generate_random_A_productState(Dvec, d);
else
    [A] = generate_random_A_productState(Dvec, d);
end

if TDVP == true
 %   mainLoop_TDVP;
    mainLoop_TDVP_Al;
else %DMRG == true
    mainLoop_DMRG;
end
    
