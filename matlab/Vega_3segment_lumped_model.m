% Vega Rocket 3 segment model
% Creates state space matrices for pitch dynamics (no aerodynamics)
% Joe T. 17/Oct/2016
%
% Coefficient matrices below for equations of the form:
%
% M*dx/dt = A*x + B*u
%       
% x = [a1, da1/dt, a2, da2/dt, a3, da3/dt]^T
% u = delta
%
% delta is gimbal angle in radians
% a1 is the rocket pitch angle
% a2 and a3 are the angles between segments 1&2 and 2&3 respectively

% Physical parameters for Vega
m = 0.5e5; % mass of segment kg
h = 10; % height of segment m
I = 5e5; % MOI of segment kg m^2
T = 2.3e6; % forward thrust N
k = 6e7; % joint stiffness Nm

% insel = 1
% insgn = 1
% outsel = 1

insel = 3
insgn = -1
outsel = 2

% insel = 4
% insgn = -1
% outsel = 3

% Assemble matrices

M1 = [   1,                 0, 0,                 0, 0,                0;
        0,                 0, 1,                 0, 0,                0;
        0,                 0, 0,                 0, 1,                0;
        0,     2*m*h^2 + 3*I, 0, (3*m*h^2)/2 + 2*I, 0,    (m*h^2)/2 + I;
        0, (3*m*h^2)/2 + 2*I, 0, (7*m*h^2)/6 + 2*I, 0, (5*m*h^2)/12 + I;
        0,     (m*h^2)/2 + I, 0,  (5*m*h^2)/12 + I, 0,    (m*h^2)/6 + I];

A1 = [0, 1,             0, 0,           0, 0;
     0, 0,             0, 1,           0, 0;
     0, 0,             0, 0,           0, 1;
     0, 0,     (2*T*h)/3, 0,     (T*h)/6, 0;
     0, 0, (2*T*h)/3 - k, 0,     (T*h)/6, 0;
     0, 0,       (T*h)/6, 0, (T*h)/6 - k, 0];
 
B1 = [  0,          0,      0,      0;
        0,          0,      0,      0;
        0,          0,      0,      0;
        (3*T*h)/2,  h/2,    -h/2,   -3*h/2;
        (2*T*h)/3,  2*h/3,  -h/3,   -4*h/3;
        (T*h)/6,    h/6,    h/6,    -5*h/6];

% State space matrices

A = -M1\A1;
B = insgn.*M1\B1;
C = [1,0,0,0,0,0;
     1,0,1,0,0,0;
     1,0,1,0,1,0];
D = 0;

R = -M1(4:6,[2,4,6])\A1(4:6,[1,3,5]);
Q = M1(4:6,[2,4,6])\B1(4:6,:);
S = C(:,[1,3,5]);

[V La] = eig(R);
Ql = inv(V)*Q
Sl = S*V

Ql = insgn.*Ql(:,insel);
Sl = Sl(outsel,:);

sca = diag(sqrt(Sl./Ql.'));
Qla = sca*Ql;
x1 = (Qla)./norm(Qla);

[J U] = lanczos(La,x1);

n=3;
b_diag = -diag(J,1);
a_diag = diag(J,0);
x = zeros(n,1);
x(1) = 1/norm(Qla);
x(2) = a_diag(1)*x(1)./b_diag(1);
x(3) = a_diag(2)*x(2)./b_diag(2)-b_diag(1)*x(1)./b_diag(2);

Di = diag(x);
Kn = Di*J*Di
Mn = Di^2
Anew = Mn\Kn;
[o3 o4] = eig(Anew);
mrat = diag(Mn)./Mn(1,1);
krat = diag(Kn,1)./Kn(1,2);

%% tests
% P1 = V*inv(sca)*U.'*Di;
% inv(P1)*R*P1;
% inv(P1)*Q;
% S*P1;
% 
% % WBC Test
% w = sym('omega');
% k0 = sym('k0');
% 
% % WTF model
% Ag = [0,1;-w^2,-w]; Bg = [0;w^2]; Cg = [1,0]; Dg = 0;
% ordg = size(Ag,1);
% 
% % WBC dd model
% Aw = [Ag,-Bg*Cg;-Bg*Cg,Ag]; Bw = [zeros(2,1),Bg;Bg,zeros(2,1)]; Cw = [Cg,zeros(1,2)];
% % WBC fd model
% P = [1/k0,1,0;0,1,0;0,0,1];
% Ax = Aw; Bx=[Bw,zeros(2*ordg,1)]*P; Cx=k0*Cw; Dx=[0,-k0,k0/2]*P;
% % WBC ideal actuator
% Bx1=Bx(:,1); Bx2=Bx(:,2:3);
% Dx1=Dx(:,1); Dx2=Dx(:,2:3);
% Ay=Ax+Bx1*Cx; By=Bx2+Bx1*Dx2; Cy=Cx; Dy=Dx2;

%%
sys_ol = ss(A,-B(:,insel),C(outsel,:),D);
s = tf('s');
wtfn = @(w,n) 2/(2 + (2*s*n/w) + (s^2*n^2/w^2));

a = 1;
n = 3;

k0 = Kn(1,1);
w = sqrt(k0/Mn(1,1));

G = wtfn(w,1);
G2np1 = wtfn(w,n);
C1 = 0.5*(1-a*G^2+(1-a)*G2np1);
C2 = a*G;

% system with added imaginary spring
sys_k0 = series(k0,feedback(sys_ol,k0));

% closed loop wbc control system
sys_cl = series(C1,feedback(sys_k0,C2,+1));

[Y,T,X] = step(sys_cl,5);
a1 = X(:,1)
a2 = X(:,1) + X(:,3)
a3 = X(:,1) + X(:,3) + X(:,5)

figure(1);
title('Step Response')
plot(T,a1,'b',T,a2,'r-.',T,a3,'k--')
axis([0 5 -0.3 1.1])
axis([0 5 -0.7 2])
legend('\alpha_1','\alpha_2','\alpha_3');
title('Input f_2, output \alpha_2')
xlabel('t')
ylabel('Segment Attitude Angles')

oldFolder = cd('K:\Documents\eccomas2017\graphics')
matlab2tikz('graphf2.tex', 'height', '5cm' , 'width', '15cm');
cd(oldFolder)