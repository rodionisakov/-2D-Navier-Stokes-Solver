function [A] = A(dx,dy,N,M)
%Function for Poisson solver - A matrix

NM = N-1*M; %Not calcuating inlet pressure so 1 less column
A = sparse(NM,NM); %sparse function reduces the computing time by squeezing out zero values 

%% Centre Points - (T1)
for i=2:N-2
    for j=2:M-1
    q=((j-1)*(N-1))+i; %pointer A
    
A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q+1) = 1/dx^2; %east coefficent
A(q,q-1) = 1/dx^2; %west coefficent
A(q,q+(N-1)) = 1/dy^2; %top coefficent
A(q,q-(N-1)) = 1/dy^2; %bottom coefficent
    end
end

%% Left Points - (T2)
i=1;
for j=2:M-1

    q=((j-1)*(N-1))+i; %pointer A

A(q,q+1) = 1/dx^2; %east coefficent
A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q+(N-1)) = 1/dy^2; %top coefficent
A(q,q-(N-1)) = 1/dy^2; %bottom coefficent

end

%% Right Points - (T3)
i=N-1;
for j=2:M-1

    q=((j-1)*(N-1))+i;%pointer A

A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q-1) = 2/dx^2; %west coefficent
A(q,q+(N-1)) = 1/dy^2; %top coefficent
A(q,q-(N-1)) = 1/dy^2; %bottom coefficent
  
end

%% Bottom Points - (T4)
j=1;
for i=2:N-2

    q=((j-1)*(N-1))+i;%pointer A

A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q+1) = 1/dx^2; %east coefficent
A(q,q-1) = 1/dx^2; %west coefficent
A(q,q+(N-1)) = 2/dy^2; %top coefficent
    
end

%% Top Points - (T5)
j=M;
for i=2:N-2

    q=((j-1)*(N-1))+i;%pointer A

A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q+1) = 1/dx^2; %east coefficent
A(q,q-1) = 1/dx^2; %west coefficent    
A(q,q-(N-1)) = 2/dy^2; %bottom coefficent
    
end

%% Bottom Left points- (T6)

    i=1;
    j=1;
    q=((j-1)*(N-1))+i;%pointer A

A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q+1) = 1/dx^2; %east coefficent
A(q,q+(N-1)) = 2/dy^2; %top coefficent

%% Bottom Right points- (T7)

    i=N-1;
    j=1;
    q=((j-1)*(N-1))+i;%pointer A

A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q-1) = 2/dx^2; %west coefficent  
A(q,q+(N-1)) = 2/dy^2; %top coefficent

%% Top Left points- (T8)

    i=1;
    j=M;
    q=((j-1)*(N-1))+i;%pointer A

A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q+1) = 1/dx^2; %east coefficent
A(q,q-(N-1)) = 2/dy^2; %bottom coefficent


%% Top Right points- (T9)

    i=N-1;
    j=M;
    q=((j-1)*(N-1))+i;%pointer A
  
A(q,q) = (-2/dy^2)-(2/dx^2); %point coefficent
A(q,q-1) = 2/dx^2; %west coefficent  
A(q,q-(N-1)) = 2/dy^2; %bottom coefficent


end