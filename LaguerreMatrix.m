function LaguerreMatrix=LaguerreMatrix(n_max)
LaguerreMatrix=zeros(n_max+1,n_max+1,n_max+1);
for k=0:n_max
    for j=0:k
    LaguerreMatrix(j+1,k+1,(j+1):-1:1)=LaguerreGen(j,k-j);
    end;
end;


