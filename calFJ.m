function [fout,jout] = calFJ(u)
global N

%two independent variables x and y:
n=N/2; %splitting two sets
x=ones(n,1); %column vector x
y=ones(n,1); %column vector y

%assigning to the input vectors of the function:
o=1:2:N; %odd
e=2:2:N; %even
x=u(o,:); %spreading x among the odd rows of u
y=u(e,:); %spreading y among the even rows of u

%Defining the functions to differentiate:
fout=ones(N,1);

%for first two rows: 
i=1;
fout(i,1)=1+x(i)^2 *y(i)-4*x(i)+(-2*x(i)+x(i+1));
fout(i+1,1)=3*x(i)-x(i)^2 *y(i)+(-2*y(i)+y(i+1));

%for last two rows:
i=N-1;
fout(i,1)=1+x(n)^2 *y(n)-4*x(n)+(x(n-1)-2*x(n));
fout(i+1,1)=3*x(n)-x(n)^2 *y(n)+(y(n-1)-2*y(n));

%for rest of the rows:
j=1;
if n>2
     for i=2:n-1
         if j<N
             j=j+2;
             fout(j,1)=1+x(i)^2 *y(i)-4*x(i)+(x(i-1)-2*x(i)+x(i+1));
             fout(j+1,1)=3*x(i)-x(i)^2 *y(i)+(y(i-1)-2*y(i)+y(i+1)); 
         end
     end
end
%Analytical Jacobian:
jout=zeros(N,N); %blank jacobian

%first two rows:
i=1;
%first row:
jout(i,i)=2*x(i).*y(i)-4-2;
jout(i,i+1)=x(i).^2;
jout(i,i+2)=1;

%second row:
j=2;
jout(j,i)=3-2*x(i).*y(i);
jout(j,i+1)=-x(i).^2-2;
jout(j,i+3)=1;

%rest of the rows:
%odd rows nonzero points:
%dfi/du(i-1):
j=0;
for i=3:2:N
    if j<N-2
        j=j+1;
        jout(i,j)=1;
        j=j+1;
    end
end
%dfi/dui:
j=1;
for i=3:2:N 
   if j<n
       j=j+1;
   jout(i,i)=2*x(j).*y(j)-4-2;
   end
end
%dfi/du(i+1):
j=3;
if N>4
   for i=3:2:N-2
       if j<N
           j=j+2;
           jout(i,j)=1;
       end
   end
end
%dfi/dvi:
j=1;
for i=3:2:N
    if j<N
        j=j+1;
        jout(i,2*j)=x(j).^2;
    end
end

%Even rows nonzero points:
%dgi/dui:
j=1;
for i=4:2:N
   if j<N
      j=j+2;
      jout(i,j)=3-2*x(i/2).*y(i/2);
   end
end
%dgi/dv(i-1):
j=0;
for i=4:2:N
   if j<N-2
       j=j+2;
   jout(i,j)=1;
   end
end
%dgi/dvi:
for i=4:2:N
    jout(i,i)=(-x(i/2).^2)-2;
end
%dgi/dv(i+1):
j=4;
if N>4
    for i=4:2:N-2
        if j<N
           j=j+2;
           jout(i,j)=1;
        end
    end
end
end
