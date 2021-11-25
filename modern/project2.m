%% Given

A = [-6.30908 0 -10.9544 0 83.74458 0 0 24.05866;
     0 -161.083 0 0 51.52923 0 -18.0261 0;
     -18.7858 0 -46.3136 0 275.6592 0 0 158.3741;
     0 0 0 -17.3506 193.9373 0 0 0;
     1.299576 0 2.969317 0.3977 -38.7024 0.105748 0 0;
     16.64244 0 38.02522 5.066579 -479.384 0 0 0;
     0 -450.386 0 0 142.2084 0 -80.9472 0;
     2.02257 0 4.621237 0 0 0 0 -51.2108];

BT = [0 0 0 3.946683 0 0 0 0];

B = BT';

C = [0 0 0 5.066579 -116.446 0 0 0; 
    0 0 0 0 1 0 0 0;
    12.96989 10.32532 -0.56926 0 0 0 0 0];

n = rank(A); % system order
m = rank(B); % # of inputs
c = rank(C); % rank of output
l = 3; % # of outputs

D = zeros(c, m);


%% Full Order Observer

x0 = [1 2 3 4 5 6 7 8]';
y0 = C*x0;
% Least-Squares Estimation of x0hat
x0hat = pinv(C'*C)*C'*y0;

dpoles = eig(A);
obspoles = linspace(-1000,-300,8);
KT = place(A', C', obspoles);
F = place(A, B, dpoles);
K = KT';

Aobs = A - K*C;
Bobs = [B K];
Cobs = eye(n);
Dobs = zeros(n, m+c);

w = eig(Aobs); % Negative
u = eig(A - B*F); % Negative as well

%% Full Order Observer Plots

% figure(1);
% plot(out.x);
% title("x from Full Order Observer");
% ylabel("x");
% 
% figure(2);
% plot(out.xhat);
% title("xhat from Full Order Observer");
% ylabel("xhat");
% 
% figure(3);
% plot(out.inerror);
% title("Input Error from Full Order Observer");
% ylabel("x - xhat");
%
% figure(4);
% plot(out.outerror);
% title("Output Error from Full Order Observer");
% ylabel("y - yhat");
% 
% figure(5);
% plot(out.y);
% title("Output Response from Full Order Observer");
% ylabel("y");

%% Reduced-Order Observer

r = n - c;
% C*L = Ic, C1*L1=Ir, C*L1 = 0, C1*L=0 -> Lest = C'*C
Lest = C'*C;
C1 = null(Lest)';
Caug = [C;C1];
Laug = inv(Caug);
L = Laug(1:n,1:c);
L1 = Laug(1:n,c+1:n);
O = obsv(C1*A*L1, C*A*L1);
% observable, rank is r = 5
disp(['Rank of Reduced-Order Observer:', num2str(rank(O))]); 

lambda_sys = eig(A - B*F);
lambda_red_obs = lambda_sys(1:5); %linspace(-1000,-300,5);
K1T = place((C1*A*L1)', (C*A*L1)', lambda_red_obs);
K1 = K1T';
Aq = C1*A*L1 - K1*C*A*L1;
Bq = C1*B - K1*C*B;
Kq = C1*A*L1*K1 + C1*A*L - K1*C*A*L - K1*C*A*L1*K1;
% inv(L1'*L1) = eye(5)
x0hat_reduced = eye(5)*L1'*(pinv(C'*C)*C'-(L+L1*K1))*y0;

%% Reduced Order Plots

figure(6);
plot(out.x_red);
title("x from Reduced Observer");
ylabel("x");

figure(7);
plot(out.xhat_red);
title("xhat from Reduced Observer");
ylabel("xhat");

figure(8);
plot(out.inerror_red);
title("Input Error from Reduced Observer");
ylabel("x - xhat");

figure(9);
plot(out.y_red);
title("Output Response from Reduced Order Observer");
ylabel("y");