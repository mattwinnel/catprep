function [ psi_out ] = scheme_II_function(n,dim,eta_total,squeezing_vacuum,squeezing_Fock,ROUND_MAX)

%
% this function numerically performs scheme II for the preparation of
% optical cat states from ''Deterministic preparation of optical cat and
% Gottesman-Kitaev-Preskill (GKP) states'' by Winnel et al. (2023).

FIRST_ANGLE = [];
% inputs
% n: size of the initial Fock state at the input
% dim: max cut_off photon number (i.e., for qubits dimension = 2, dim = 1)
% eta_total: total transmissivity = eta^k
% squeezing: inline squeezing in dB
% ROUND_MAX: maximum number of rounds = k

% psi_out: pure output state

dB = squeezing_vacuum;
r_vacuum = dB/(20*log10(exp(1))); % r is the squeezing parameter with variance exp(2r);

dB = squeezing_Fock;
r_Fock = dB/(20*log10(exp(1))); % r is the squeezing parameter with variance exp(2r);


a = circshift(diag(sqrt(0:1:dim)),-1); % annihilation operator


I = eye(dim+1); % identity matrix


eta = eta_total^(1/ROUND_MAX); % transmissivity eta such that total transmissivity = eta^k;


% beamsplitter Gaussian unitary
theta = acos(sqrt(eta));
B = expm(theta*(kron(a',a)-kron(a,a')));




% choose a random angle from the circle discritised into k angles
A = linspace(pi,0,ROUND_MAX+1);
A = A(2:end);
III = randperm(length(A));
A = A(III);


% initially the input state (called psi_out) is a Fock state |n>
% psi_out will evolve to eventually be the output state
psi_out = zeros(1,dim+1).';
psi_out(n+1) = 1;
% squeeze initial Fock state
% squeezing Gaussian unitary
S_Fock = expm(0.5*(r_Fock+0.00000000000001)*(a^2-(a')^2));
psi_out = S_Fock*psi_out;


% define vacuum state
fock_0 = zeros(dim+1,1);
fock_0(1) = 1;



for ROUND = 1:ROUND_MAX%:10
    
    % inititial state to the ROUND is equal to the final state of
    % the previous round
    psi = psi_out;
    
    % angle of squeezing
    theta = A(ROUND);
  
    
    
    % squeeze vacuum
    % squeeze mode to be measured at angle theta:
    % squeezing Gaussian unitary
    S_vacuum = expm(0.5*r_vacuum*(a^2*exp(-2*1i*theta)-(a')^2*exp(2*1i*theta)));    
    

    % interact the input mode with vacuum on a beamsplitter
    psi_out = B*kron(S_vacuum*fock_0,psi);
    
    
    

    
    
    % randomly select the measurement outcome from the probability
    % distribution
    rho_trace = TrX(psi_out*psi_out',[2],[dim+1 dim+1]);
    R = 0;% Don't squeeze measurement
    [ m ] = histogram_plots( rho_trace/trace(rho_trace),dim,R );

    
    
    % measurement operator
    fock_m = zeros(dim+1,1);
    fock_m(m+1) = 1;
    M = (fock_m); % don't squeeze measurement
    

    % perform measurement
    psi_out = kron(M,I)'*psi_out;
    

    % renormalise state
    psi_out = psi_out/norm(psi_out);

    
end

end

