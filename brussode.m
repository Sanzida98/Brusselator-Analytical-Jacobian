
% if nargin<1 %number of function input
%    N = 10; %grid points
% end
N=100;
tspan = [0; 10]; %time interval
y0 = [1+sin((2*pi/(N+1))*(1:N)); repmat(3,1,N)]; %initial condition: u(j)=1+sin(2*pi*x) and v(j)=0
options = odeset('Vectorized','on','JPattern',jpattern(N)); %vectorizing function output and taking jpattern from jpattern(N)
[t,y] = ode15s(@f,tspan,y0,options); %ode15s call
y;
u = y(:,1:2:end);
x = (1:N)/(N+1); %defining x
figure; 
surf(x,t,u); %3D x:space, t:time
view(-40,50);
xlabel('space');
ylabel('time');
zlabel('solution u');
title(['The Brusselator for N = ' num2str(N)]);
save('C:\Users\saeid\Desktop\SHORNA\UMBC_CPL\Brusselator\Jacobian\jverify\jverify\variable.mat','y');

%function declare:
    function dydt = f(t,y) %function declare
      % Derivative function
      N = 100;
      c = 0.02 * (N+1)^2; 
      dydt = zeros(2*N,size(y,2));      % preallocate dy/dt
      
      % Evaluate the 2 components of the function at one edge of the grid
      % (with edge conditions).
      i = 1;
      dydt(i,:) = 1 + y(i+1,:).*y(i,:).^2 - 4*y(i,:) + c*(1-2*y(i,:)+y(i+2,:));
      dydt(i+1,:) = 3*y(i,:) - y(i+1,:).*y(i,:).^2 + c*(3-2*y(i+1,:)+y(i+3,:));
      % Evaluate the 2 components of the function at all interior grid points.
      i = 3:2:2*N-3;
      dydt(i,:) = 1 + y(i+1,:).*y(i,:).^2 - 4*y(i,:) + c*(y(i-2,:)-2*y(i,:)+y(i+2,:));
      dydt(i+1,:) = 3*y(i,:) - y(i+1,:).*y(i,:).^2 + c*(y(i-1,:)-2*y(i+1,:)+y(i+3,:));
      
      % Evaluate the 2 components of the function at the other edge of the
      % grid
      i = 2*N-1;
      dydt(i,:) = 1 + y(i+1,:).*y(i,:).^2 - 4*y(i,:) + c*(y(i-2,:)-2*y(i,:)+1);
      dydt(i+1,:) = 3*y(i,:) - y(i+1,:).*y(i,:).^2 + c*(y(i-1,:)-2*y(i+1,:)+3);
    end
  % brussode
%Subfunction -- the sparsity pattern
function S = jpattern(N)
% Jacobian sparsity pattern
N=100;
B = ones(2*N,5);
B(2:2:2*N,2) = zeros(N,1);
B(1:2:2*N-1,4) = zeros(N,1);
S = spdiags(B,-2:2,2*N,2*N);
end
