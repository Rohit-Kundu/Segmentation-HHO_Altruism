% Main paper:
% Segmentation of Brain MRI using an Altruistic Harris Hawks' Optimization algorithm
% Rajarshi Bandopadhyay, Rohit Kundu, Diego Oliva, Ram Sarkar
% _____________________________________________________

function [X]=initialization(N,dim,up,down)

if size(up,1)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,1)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
    %X(1,:)=;random one
    for i = 2:N
        %sine
        %X(i,:)=sin(pi*X(i-1,:));
        
        %singer
        %X(i,:)=1.07*(7.868*(X(i-1,:)) -23.31*(X(i-1,:)*X(i-1,:))+ 28.75*(X(i-1,:)*X(i-1,:)*X(i-1,:))- 13.302875*(X(i-1,:)*X(i-1,:)*X(i-1,:)*X(i-1,:)));
        
        %sinusoidal
        %X(i,:)= 2.3*(X(i-1,:)*X(i-1,:))* sin(pi*X(i-1));
        
        %chebyshev
        %X(i,:)=cos(acosd(X(i-1,:)));
        
        %Tent
        %if X(i-1,:)<0.7
        %    X(i,:)= X(i-1,:)/0.7;
        %end
        %if X(i-1,:)>=0.7
        %    X(i,:)= (10/3)*(1-X(i-1,:));
        %end
        
        %Logistic
        X(i,:)= 4*X(i-1,:)*(1-X(i-1,:));
        
        %Iterative
        %X(i,:)=sin(0.7*pi/(X(i-1,:)));
        
        %Gauss
        %if X(i-1,:)==0
        %    X(i,:)= 1
        %end
        %if X(i-1,:)~=0
        %   X(i,:)= 1/mod(X(i-1,:))
        %end
    end
end
