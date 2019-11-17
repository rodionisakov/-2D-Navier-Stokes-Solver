function [U] = U(W1,c,k,P,dx,dy,N,M,VOF)
%Correction U velocity

U = W1; %so that the outlet is not calculated

%% i=1 Inlet

    U(1,:) = 1;
    
%% Whole body inc top and bottom

for i=2:N-1
    for j=1:M
        
U(i,j) = W1(i,j)-(1/c(k))*((P(i+1,j)-P(i-1,j))/(2*dx)); %U point value

    end
end

U = U.*VOF; %remove velocity from inside square cylinder

end

