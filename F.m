function [F1C,F2C,F1V,F2V] = F(F1C,F2C,F1V,F2V,Re,u,v,dx,dy,N,M)
%Function Convective and Velocity flux

%% Inlet boundary condition
i=1;
for j = 2:M-1
F1C(i,j) = -((u(i+1,j)^2 - (2-u(i+1,j))^2)/(2*dx)) - (((u(1,j)*v(i,j+1))-(u(i,j-1)*v(i,j-1)))/(2 * dy));%convective flux
F2C(i,j) = -((u(i+1,j)*v(i+1,j))-(-(2-u(i+1,j)*v(i+1,j)))/(2 * dx)) - ((v(i,j+1)^2 - v(i,j-1)^2)/(2 * dy));%convective flux
F1V(i,j) = (1/Re)*(((u(i+1,j)-u(i+1,j))/dx^2)+((u(i,j+1)+1+u(i,j-1))/dy^2));%viscous flux
F2V(i,j) = (1/Re)*((v(i,j+1)+v(i,j-1))/dy^2);%viscous flux

end
%% Centre points
for i = 2:N-1
    for j = 2:M-1  
F1C(i,j) = -((u(i+1,j)^2 - u(i-1,j)^2)/(2 * dx)) - (((u(i,j+1)*v(i,j+1))-(u(i,j-1)*v(i,j-1)))/(2 * dy));%convective flux
F2C(i,j) = -(((u(i+1,j)*v(i+1,j))-(u(i-1,j)*v(i-1,j)))/(2 * dx)) - ((v(i,j+1)^2 - v(i,j-1)^2)/(2 * dy));%convective flux
F1V(i,j) = (1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2));%viscous flux
F2V(i,j) = (1/Re)*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2)+((v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2));%viscous flux
    end
end
%% Top boundary
 j = M;
for i = 2:N-1 
F1C(i,j) = -((u(i+1,j)^2 - u(i-1,j)^2)/(2 * dx)) + ((u(i,j-1) * v(i,j-1))/dy);%convective flux
F1V(i,j) = (1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+(2*(u(i,j-1)-u(i,j))/dy^2)); %viscous flux
end
%% Bottom boundary condition
j = 1; 
for i = 2:N-1
F1C(i,j) = -((u(i+1,1)^2 - u(i-1,1)^2)/(2 * dx))+ ((u(i,j+1) * v(i,j+1))/dy);%convective flux
F1V(i,j) = (1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+(2*(u(i,j+1)-u(i,j))/dy^2));%viscous flux
end
%% Outlet boundary condition
%orlansky boundary condition viscous flux is absent
i=N;
for j=2:M-1
F1C(i,j) = -((3*u(i,j)-4*u(i-1,j)+u(i-2,j))/(2*dx)) - (((u(i,j+1)*v(i,j+1))-(u(i,j-1)*v(i,j-1)))/(2 * dy));%convective flux
F2C(i,j) = -((3*u(i,j)*v(i,j)-4*u(i-1,j)*v(i-1,j)+u(i-2,j)*v(i-2,j))/(2*dx))- ((v(i,j+1)^2 - v(i,j-1)^2)/(2 * dy));%convective flux
end 
%% Top right point 
%orlansky boundary condition, viscous flux is absent
F1C(N,M) = -(3*u(N,M)-4*u(N-1,M)+u(N-2,M))/(2*dx) + (u(N,M-1) * v(N,M-1))/dy;%convective flux
F2C(N,M) = -((3*u(i,j)*v(i,j)-4*u(i-1,j)*v(i-1,j)+u(i-2,j)*v(i-2,j))/(2*dx));%convective flux
%% Bottom right point
%orlansky boundary condition, viscous flux is absent
F1C(N,1) = -(3*u(N,1)-4*u(N-1,1)+u(N-2,1))/(2*dx) - (u(N,1+1) * v(N,1+1))/dy;%convective flux
F2C(N,1) = -((3*u(i,j)*v(i,j)-4*u(i-1,j)*v(i-1,j)+u(i-2,j)*v(i-2,j))/(2*dx));%convective flux

end