function [ m,peak,y_max ] = histogram_plots( rho,dim,r)
X=0:dim-1;
p=[];

a=(circshift(diag(sqrt(0:1:dim)),-1));
S = expm(0.5*r*(a^2-(a')^2));
for m = X
    fock_m = zeros(dim+1,1);
    fock_m(m+1) = 1;
    fock_m = S*fock_m;
    M = fock_m*fock_m';
    p = [p trace(M*rho)];
end

VECTOR = zeros;
for j = 1:length(X)
    N = round(p(j)*100000);
    VECTOR = [VECTOR ones(1,N).*X(j)];
end
%VECTOR
m = randsample(VECTOR,1);


[y_max, index] = max(p);
X_max = X(index);
peak = X_max;



end

