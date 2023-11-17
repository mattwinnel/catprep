
dim = 30; % dim is the truncation Fock number, i.e., dim = 1 for qubits

% scheme I

n = 14;
ROUND_MAX = 6;
squeezing = 8; % dB.
eta_total = 0.25;


[ psi_out_scheme_I ] = scheme_I_function(n,dim,eta_total,squeezing,ROUND_MAX);



dim_trun = 16;
figure
wigner_function2D(psi_out_scheme_I(1:dim_trun+1)*psi_out_scheme_I(1:dim_trun+1)',dim_trun)
drawnow
disp('displaying output cat for scheme I')
disp('press any key to continue')
pause

n = 14;
ROUND_MAX = 100;
squeezing_vacuum = 6; % dB.
squeezing_Fock = 0; % dB.
eta_total = 0.5;



[ psi_out_scheme_II ] = scheme_II_function(n,dim,eta_total,squeezing_vacuum,squeezing_Fock,ROUND_MAX);

figure
wigner_function2D(psi_out_scheme_II(1:dim_trun+1)*psi_out_scheme_II(1:dim_trun+1)',dim_trun)
drawnow
disp('displaying output cat for scheme II')



