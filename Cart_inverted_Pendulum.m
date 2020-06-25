%% Defining Parameters for Cart inverted pendulum
clear;
close all;
s=tf('s');
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
% M = .5;
% m = 0.2;
% b = 0.1;
% I = 0.006;
% g = 9.8;
% l = 0.3;
% 
% p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices
% 
% A = [0      1              0           0;
%      0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
%      0      0              0           1;
%      0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
% B = [     0;
%      (I+m*l^2)/p;
%           0;
%         m*l/p];
% C = [1 0 0 0;
%      0 0 1 0];
% D = [0;
%      0];

% http://databookuw.com/databook.pdf
m = 1; M = 5; L = 2; g = -10; d = 1;
b = 1; % Pendulum up (b=1)
A = [0 1 0 0;
0 -d/M b*m*g/M 0;
0 0 0 1;
0 -b*d/(M*L) -b*(m+M)*g/(M*L) 0];
B = [0; 1/M; 0; b*1/(M*L)];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,...
'outputname',outputs);

%% Initial stability analysis
% To check location of poles
poles = eig(sys_ss);

% Controllability
ctrl_mat = [B A*B (A^2)*B (A^3)*B];
rank(ctrl_mat);

% Observability
obsv_mat = [C; C*A; C*(A^2); C*(A^3)];
rank(obsv_mat);

% TF Matrix
tf_matrix = C*inv(s*eye(4)-A)*B + D;
tf_matrix_min = minreal(tf_matrix);
w=logspace(-3,3,1000);

%% Controller design - LQR
R = 0.1; %Sensor variance
Q = [30 0 0 0;0 70 0 0;0 0 10 0;0 0 0 40]; % State penalty
% Q = C'*C;
G = lqr(A,B,Q,R);
Glqr=ss(A,B,G,0);

So1=inv(eye(size(Glqr))+Glqr);
To1=eye(size(Glqr))-So1;%Te=PK*inv[I+PK]

%figure;
% sv=sigma(So1,w);
% tsv=20*log10(sv);
% semilogx(w,tsv);
% title('Sensitivity (S) using LQR');
% grid on;
% xlabel('Frequency (rad/sec)');
% ylabel('Singular Values (dB)');

%figure;
% sv=sigma(To1,w);
% tsv=20*log10(sv);
% semilogx(w,tsv);
% title('Complementary Sensitivity (T) using LQR');
% grid on;
% xlabel('Frequency (rad/sec)');
% ylabel('Singular Values (dB)');

figure;
set(gcf,'Position',[1000 350 500 400])
step(To1);
title('Step response using LQR without prefilter');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Magnitude');
hold on;
stepinfo(To1)

%% Simulate cart using LQR controller
tspan = 0:.1:20;
x0 = [-1; 0; pi-0.5; 0]; % initial condition
wr = [1; 0; pi+0; 0]; % reference position
u=@(x)-G*(x - wr); % control law
[t,y] = ode45(@(t,x)cartpend(x,m,M,L,g,d,u(x)),tspan,x0);

figure;
for k=1:length(t)
    drawcartpend(y(k,:),m,M,L);
end
