s = tf('s');
t = 0:0.001:1;

%% Question A %%

% Defining transfer function G(s) and plotting root locus
k1 = 1;
z1 = [inf;inf;inf];
p1 = [0;-10;-50];
[n1,d1] = zp2tf(z1,p1,k1);
G = tf(n1,d1);

figure;
rlocus(G);
title("Open Loop Root Locus of G(s)")

% Parameters
mpos = 0.075;
ts = 0.4;

zeta = sqrt(log(mpos)^2 / (pi^2 + log(mpos)^2));
natfreq = 3/(zeta * ts);
x = natfreq*zeta; % b characteristic
y = natfreq*sqrt(1 - zeta^2); % Damping frequency

% Poles
lambda1 = x + y*1i;
lambda2 = x - y*1i;

% Find controller zero zA
theta50 = atand(y/(50 - x));
theta10 = atand(y/(10 - x));
theta0 = 180 - atand(y/x);
thetap = theta50 + theta10 + theta0;
thetaA = -1*180 + thetap;
zA = y/tand(thetaA) + x;

% Find controller gain K
L50 = sqrt((50 - x)^2 + y^2);
L10 = sqrt((10 - x)^2 + y^2);
L0 = sqrt((x)^2 + y^2);
LA = sqrt((zA - x)^2 + y^2);
K = (L50*L10*L0)/LA;

% PD contoller
numPD = [ K K*zA ];
demPD = 1;
Gpd = tf(numPD, demPD);

% Open Loop
Gcpd = series(Gpd,G);
[nG, dG] = tfdata(Gcpd,'v');
figure;
rlocus(Gcpd);
title("Open Loop PD Controller Root Locus")

% Closed Loop
[nGf, dGf] = feedback(nG, dG, 1, 1, -1);
Gf = tf(nGf,dGf);
pG = pole(Gf);
zG = zero(Gf);
figure;
rlocus(Gf);
title("Closed Loop PD Controller Root Locus")

% Step Response
figure;
step(Gf, t);
title('Step Response of PD Controller')

%% Question B %%

% Ramp Steady State Error
% Since G(s) is a type1 system, ramp steady state equals 1/Kv
Kv = 1/(50*10);
ess = 1/Kv; % evaluates to 500

% PI Controller
numPI = [ 1 0.01 ];
demPI = [ 1 0 ];
Gpi = tf(numPI,demPI);
Gc = series(Gpi,Gcpd);
[nGc,dGc] = tfdata(Gc,'v');
[nGclose, dGclose] = feedback(nGc, dGc, 1, 1, -1);
Gclosed = tf(nGclose, dGclose);
figure;
rlocus(Gclosed);
title('Controller R-locus')

% Step Response
figure;
step(Gclosed,t);
title('PID Controller Step Response');

% Ramp Response
figure;
plot(t,t,'--');
hold on;
step(Gclosed/s, t);
title('PID Controller Ramp Response');

%% Question C %%

% Phase-Lead Controller
% Since the available contributed phase is less than zero (thetac = -36.31), we can cancel
% out one of the non-dominant poles, so the phase goes up to a positive
% value for the pole to contribute to. Arbitrarily, if zl = 10, then pl is
% below:

zlead = 10;
thetapl = 180 - (theta50 + theta0);
px = y/tand(thetapl) + x;
Lpl = sqrt((px - x)^2 + y^2);
K2 = L50*L0*Lpl;

nplead = [K2 K2*zlead];
dplead = [1 px];
plead = tf(nplead,dplead);

% Open Loop
Gopen = series(plead,G);
[nGopen, dGopen] = tfdata(Gopen,'v');
figure('Name', 'Open Loop Lead R-locus');
rlocus(Gopen);
title('Open Loop Lead R-locus');

% Closed Loop
[nLead, dLead] = feedback(nGopen,dGopen,1,1,-1);
Glead = tf(nLead,dLead);
pGlead = pole(Glead);
zGlead = zero(Glead);
figure('Name', 'Closed Loop Lead R-locus');
rlocus(Glead);
title('Closed Loop Lead R-locus')

% Step Response
figure('Name', 'Lead Step Response');
step(Gf, t);
hold on;
step(Glead, t);
title('PD and Lead Controller Closed-Loop Step Responses');

% Ramp Response
figure('Name', 'Lead Ramp Response');
plot(t,t,'--');
hold on;
step(Glead/s, t);
title('Lead Ramp Response');

%% Question D %%
% Given a type 1 system, the error is proportional to the desired ramp steady
% state error by 1/Kv
syms N(w) w
[limn, limd] = tfdata(s*Gopen,'v');
N(w) = simplify(poly2sym(limn,w)/poly2sym(limd,w));
Kv = vpa(N(0));
ess = double(1/Kv);

% once we add the controller below, ess is decreased by a factor of
% 1/(1000*ess) giving an improved steady state error of 0.001. We could
% decrease pc even close to 0, but avoid zero so system type doesn't change.
nlag = [1 ess];
dlag = [1 0.001];
Glag = tf(nlag,dlag);

% Open Loop
Gc = series(Glag,Gopen); % connects lead, lag, and G
[nGc, dGc] = tfdata(Gc,'v');

% Closed Loop
[nLL, dLL] = feedback(nGc,dGc,1,1,-1);
GLL = tf(nLL,dLL);
pGLL = pole(GLL);
zGLL = zero(GLL);
figure('Name', 'Closed-loop Lead-lag R-locus');
rlocus(GLL);
title('Closed-loop Lead-lag R-locus');

% Step Response
figure('Name', 'Leag-Lag Step Response');
step(GLL, t);
title('Closed-Loop Leag-Lag Step Response');

% Ramp Response
figure('Name', 'Leag-Lag Ramp Response');
plot(t,t,'--');
hold on;
step(GLL/s, t);
title('Closed-Loop Leag-Lag Ramp Response');

