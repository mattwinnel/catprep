
dim = 30; % dim is the truncation Fock number, i.e., dim = 1 for qubits

% iterative scheme I

n = 14; % input Fock state size
ROUND_MAX = 6; % number of rounds k
squeezing = 8; % amount of squeezing in dB
eta_total = 0.25; % total transmissivity experienced by the input Fock state

% compute the output state vector for iterative scheme I
[ psi_out_scheme_I ] = scheme_I_function(n,dim,eta_total,squeezing,ROUND_MAX);


% plot Wigner function of output state
dim_trun = 16; % truncate the output state for faster Wigner function calculation
figure
wigner_function2D(psi_out_scheme_I(1:dim_trun+1)*psi_out_scheme_I(1:dim_trun+1)',dim_trun)
drawnow
disp('displaying output cat for scheme I')
disp('press any key to continue')
pause


% iterative scheme II


n = 14; % input Fock state size
ROUND_MAX = 100; % number of rounds k
squeezing_vacuum = 6; % amount of squeezing for squeezed vacuum states in dB
squeezing_Fock = 0; % amount of initial squeezing on the input state in dB
eta_total = 0.5; % total transmissivity experienced by the input state


% compute the output state vector for iterative scheme II
[ psi_out_scheme_II ] = scheme_II_function(n,dim,eta_total,squeezing_vacuum,squeezing_Fock,ROUND_MAX);

% note that scheme II works for any input state with odd/even Fock-number parity - make changes to "scheme_II_function" to try

% plot Wigner function of output state
% truncate the output state for faster Wigner function calculation
figure
wigner_function2D(psi_out_scheme_II(1:dim_trun+1)*psi_out_scheme_II(1:dim_trun+1)',dim_trun)
drawnow
disp('displaying output cat for scheme II')



