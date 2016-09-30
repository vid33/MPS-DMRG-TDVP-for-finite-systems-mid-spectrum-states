if ISING == true
    fprintf('Ising model, %J=%d, h=%d\n', J(2), h(2));
elseif HEISENBERG == true
    fprintf('Heisenberg model, Jx=%d, Jy=%d, Jz=%d, h=%d\n', Jx(2), Jy(2), Jz(2), h(2));
end

fprintf('At step %d, N is %d, Dmax is %d\n', currentStep, N, Dmax);
if TDVP == true
    fprintf('Norm of tangent vector is %d\n', norm_dA);
    fprintf('Energy, cacluated from F: %d\n', Fl{end});
end
fprintf('Energy, calculated directly: %d\n', energy1S + energy2S);
fprintf('Energy per site: %d\n', (1/N)*(energy1S + energy2S));
fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n');