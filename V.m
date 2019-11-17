function [V] = V(W2,c,k,P,dx,dy,N,M,VOF)
%Correction V velocity

V = W2; %so that the outlet is not calculated

%% Centre points

for i=2:N-1
    for j=2:M-1
        
V(i,j) = W2(i,j)-(1/c(k))*((P(i,j+1)-P(i,j-1))/(2*dy)); %V point value

    end
end

%% i=1 Inlet

V(1,:) = 0;

%% j=M Top

V(:,M) = 0;

%% j=1 Bottom

V(:,1) = 0;

V = V.*VOF;%remove velocity from inside square cylinder

end
