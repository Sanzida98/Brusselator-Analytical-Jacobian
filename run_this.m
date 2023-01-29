%Demonstration of Jacobian verification through second order convergence
clearvars
clc
global N
global t
T = 50;
N = 200; %no of independent vaqriables
naxis = (-N/2:1:N/2-1)';
dt = 2*T/N;
t = dt*naxis;
u=ones(N,1);
load("variable.mat");
u = y(end,:).';
% u(1,:)=1.5406;
% u(2,:)=3.0;
% u(3,:)=1.9096;
% u(4,:)=3.0;
% u(5,:)=1.9898;
% u(6,:)=3.0;
% u(7,:)=1.7557;
% u(8,:)=3.0;
% u(9,:)=1.2817;
% u(10,:)=3.0;
% u(11,:)=0.7183;
% u(12,:)=3.0;
% u(13,:)=0.2443;
% u(14,:)=3.0;
% u(15,:)=0.0102;
% u(16,:)=3.0;
% u(17,:)=0.0904;
% u(18,:)=3.0;
% u(19,:)=0.4594;
% u(20,:)=3.0;
%u = 2 + sech(t); %initial condition
flag = jverify(u);
[~,J1] = calFJ(u); %the Jacobian matrix
E=eig(J1)
figure;
scatter(real(E),imag(E))