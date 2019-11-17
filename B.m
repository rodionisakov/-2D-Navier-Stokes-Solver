function [B1] = B(PO,c,k,W1,W2,dx,dy,N,M)
%Function for Poisson solver - B matrix
%Changed Wi i index by +1 as to offset the lack of inlet pressure
%calculation

NM = N-1*M; %Not calcuating inlet pressure so 1 less column
B1 = zeros(N-1,M); %Prealocate B

%% Centre Points - (T1)
for i=2:N-2
    for j=2:M-1

B1(i,j) = c(k)*(((W1(i+2,j)-W1(i,j))/(2*dx))+((W2(i+1,j+1)-W2(i+1,j-1))/(2*dy))); %B point value

    end
end

%% Left Points - (T2)
    i=1;
for j=2:M-1

B1(i,j) = (-PO/dx^2)+c(k)*(((W1(i+2,j)-1)/(2*dx))+((W2(i+1,j+1)-W2(i+1,j-1))/(2*dy))); %B point value

end

%% Right Points - (T3)
    i=N-1;
for j=2:M-1
  
B1(i,j) = c(k)*(((W1(i-1,j)-4*W1(i,j)+3*W1(i+1,j))/(2*dx))+((W2(i+1,j+1)-W2(i+1,j-1))/(2*dy))); %B point value

end

%% Bottom Points - (T4)
    j=1;
for i=2:N-2
    
B1(i,j) = c(k)*(((W1(i+2,j)-W1(i,j))/(2*dx))+((W2(i+1,j+1))/dy)); %B point value

end

%% Top Points - (T5)
    j=M;
for i=2:N-2
    
B1(i,j) = c(k)*(((W1(i+2,j)-W1(i,j))/(2*dx))-((W2(i+1,j-1))/dy)); %B point value

end

%% Bottom Left - (T6)

    i=1;
    j=1;

B1(i,j) = (-PO/dx^2)+c(k)*(((W1(i+2,j)-1)/(2*dx))+(W2(i+1,j+1)/dy)); %B point value

%% Bottom Right - (T7)

    i=N-1;
    j=1;

B1(i,j) = c(k)*(((W1(i-1,j)-4*W1(i,j)+3*W1(i+1,j))/(2*dx))+(W2(i+1,j+1)/dy)); %B point value

%% Top Left - (T8)

    i=1;
    j=M;

B1(i,j) = (-PO/dx^2)+c(k)*(((W1(i+2,j)-1)/(2*dx))+(-W2(i+1,j-1)/dy)); %B point value


%% Top Right - (T9)

    i=N-1;
    j=M;

B1(i,j) = c(k)*((W1(i-1,j)-4*W1(i,j)+3*W1(i+1,j))/(2*dx)+(-W2(i+1,j-1)/dy)); %B point value

end