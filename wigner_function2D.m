function [W] = wigner_function2D(rho,dim)

% The maximum photon number in the density matrix

num_max=size(rho,1)-1;
% chose range %%
Range=4;
x_max=Range;
x_min=-Range;
p_max=Range;
p_min=-Range;

% choose the number of points across each axis %%
numpoint = 100;
x=x_min:((x_max-x_min)/numpoint):x_max;
p=p_min:((p_max-p_min)/numpoint):p_max;
%
% Generate the Wigner function %%

%If there is enough memory (e.g. > 1GB), vectorize K by k->K ,remove the
%for loop, and uncomment the following line:
%[X,P,M,N,K] =ndgrid(x,p,1:(num_max+1),1:(num_max+1),1:(num_max+1));

[X,P,M,N] =ndgrid(x,p,1:(num_max+1),1:(num_max+1));
X = X*sqrt(2);
P = P*sqrt(2);
L=LaguerreMatrix(num_max);
W=zeros(length(x),length(p), num_max+1, num_max+1);
A=zeros(length(x),length(p), num_max+1, num_max+1);
for k=1:num_max+1
W=W+(2-(M==N)).*(-1).^(M-1)./pi.*realsqrt(factorial(M-1)./factorial(N-1))...
    .*(2.*X.^2+2.*P.^2).^(abs(N-M)./2)...
    .*L(M+(num_max+1).*(N-1)+(num_max+1)^2.*(k-1))...
    .*(2.*X.^2+2*P.^2).^(k-1).*exp(-X.^2-P.^2)...
    .*real(rho(M+(num_max+1).*(N-1)).*exp(1i.*(N-M).*atan2(X,P)));
%display('percent loaded')
%k
%display('%')
%disp(num2str(100*k./(num_max+1)) '% loaded')
disp(['Wigner function loading: ' num2str(round(100*k./(num_max+1))) ' %']) 
end
W3D=sum(sum(W,4),3);
W=real(W3D);
W = W./2; % normalise

%To vectorize K, comment the previous and uncomment the following:
%W3D=sum(sum(sum(W,5),4),3);
%if abs(max(max(W)))>= abs(min(min(W)));
%    W(1,1) = max(max(W)); W(1,2) = -max(max(W));%-2*abs(max(max(W))-min(min(W)));
%else abs(max(max(W)))<= abs(min(min(W)));
%   W(1,1) = min(min(W)); W(1,2) = -min(min(W));%+2*abs(max(max(W))-min(min(W)));
%end

    

[P,X]=ndgrid(x,p);
%figure(103)
% surf(X,P,W3D)
% hold on
% shading interp
% colormap(redblue)
% colormap(parula(10));
% light
% lighting phong


%W3D=W;
%drawnow
X =X(1,:);
P = P(:,1);

%figure
imagesc(X,P,real(W))
caxis([-abs(max(max(W))) abs(max(max(W))) ])
abs(min(min(W)))
abs(max(max(W))) 

%imagesc([-numpoint/2 numpoint/2],[-numpoint/2 numpoint/2],real(W3D))

%X=W3D;
       % min(min(W))
      %  pause
    % Figures in thesis:
    
LINE_WIDTH = 1;
MARKER_SIZE = 8;
aspect_ratio = 1; % set to be 1 since it is a contour plot
x_length = 6.75;
x_length = 6.75*0.5;
y_length = x_length*aspect_ratio;
    FONT_SIZE_X_LABEL=12;
    FONT_SIZE_Y_LABEL=FONT_SIZE_X_LABEL;
    FONT_SIZE_Z_LABEL=FONT_SIZE_X_LABEL;
    
    
    
xlabel('$q$','Interpreter','latex','FontSize',FONT_SIZE_X_LABEL)
ylabel('$p$','Interpreter','latex','FontSize',FONT_SIZE_Y_LABEL)
%zlabel('$W(q,p)$','Interpreter','latex','FontSize',FONT_SIZE_Z_LABEL)
      
set(gcf,'color','w');
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperSize', [x_length y_length]); % aspect ratio is 3:2
set(gcf, 'PaperPosition', [0 0 x_length y_length]); % aspect ratio is 3:2
set(gca,'fontsize',FONT_SIZE_X_LABEL)
set(gca,'TickLabelInterpreter', 'latex');
colormap(redblue)

set(gca,'YDir','normal')

end

