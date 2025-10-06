%% LQR full car. I run it in continuous time and I compute the power
%% because in discrete time, there are some problems of full car's matrices discretization

clear all
clc;


yalmip('clear')

% 
ms  = 1274; %  Sprung mass  [Kg] 
muf = 35.5; % 35.5 front wheel mass [Kg]
mur = 35.5; % 35.5 rear wheel mass [Kg]
Ip  = 1523; % 1523 pitch inertia momento di inerzia attorno all'asse di becheggio [Kg*m^2]
Ir  = 606.1; % 606.1 roll inertia [Kg*m^2]
kf  = 27000; % 27000 Front suspension stiffness[N/m]
kr  = 30000; % 30000 Rear suspension stiffness [N/m]
t_f = 0.7; % distanza tra la fiancata sinistra e il centro di massa [m]
tr  = 0.75; % distanza tra la fiancata destra e il centro di massa [m]
ktf=228000; % 228000 Front tyre stiffness [N/m]
ktr=228000; % 228000 Rear tyre stiffness [N/m]
bf  = 0; % bf=1140;
br  = 0; % br=1250;
a   = 1.4; % distanza tra la punta e il centro di massa [m]
b   =1.47; % distanza tra la coda e il centro di massa [m]

% Generator's parameters

La=0.003;
Kv=48;
Ra=6;

% road signal
 n00=0.033;% 
n0=0.1;% reference spatial frequency page 5, same as the road1
Gn0=80*10^-6;% this value is for a C road
pi=180;
alpha=-2*pi*n00;
beta=2*pi*n0*sqrt(Gn0);
u=70*1000/3600;
load('W.mat');

%% Full car model parameters
set_p=[ms;muf;mur;Ip;Ir;kf;kr;t_f;tr;ktf;ktr;bf;br;a;b];

set_a=[La;Kv;Ra];

% Submatrices
A12 = [t_f a 1; -t_f a 1; tr -b 1; -tr -b 1];
A21 = [-kf*t_f/Ir kf*t_f/Ir -kr*tr/Ir -kr*tr/Ir; -kf*a/Ip -kf*a/Ip kr*b/Ip kr*b/Ip; -kf/ms -kf/ms -kr/ms -kr/ms];
A31 = diag([kf/muf,kf/muf,kr/mur,kr/mur]);
A33 = -diag([ktf/muf,ktf/muf,ktr/mur,ktr/mur]);
B21 = [t_f*Kv/Ir -t_f*Kv/Ir tr*Kv/Ir -tr*Kv/Ir; a*Kv/Ip a*Kv/Ip -b*Kv/Ip -b*Kv/Ip; Kv/ms Kv/ms Kv/ms Kv/ms];
B41 = -diag([Kv/muf,Kv/muf,Kv/mur,Kv/mur]);

% State Matrices
A = [zeros(4,4) A12 zeros(4,4) -eye(4); A21 zeros(3,3) zeros(3,4) zeros(3,4); zeros(4,4) zeros(4,3) zeros(4,4) eye(4); A31 zeros(4,3) A33 zeros(4,4)];
B = [zeros(4,4); B21; zeros(4,4);B41];
E = [zeros(4,4); zeros(3,4);-eye(4);zeros(4,4)];

% Submatrices
C11 = [-kf*t_f^2/Ir-kf*a^2/Ip-kf/ms,kf*t_f^2/Ir-kf*a^2/Ip-kf/ms,-kr*t_f*tr/Ir+kr*a*b/Ip-kr/ms,kr*t_f*tr/Ir+kr*a*b/Ip-kr/ms;...
        kf*t_f^2/Ir-kf*a^2/Ip-kf/ms,-kf*t_f^2/Ir-kf*a^2/Ip-kf/ms,kr*t_f*tr/Ir+kr*a*b/Ip-kr/ms,-kr*t_f*tr/Ir+kr*a*b/Ip-kr/ms;...
       -kf*t_f*tr/Ir+kf*a*b/Ip-kf/ms,kf*t_f*tr/Ir+kf*a*b/Ip-kf/ms,-kr*tr^2/Ir-kr*b^2/Ip-kr/ms,kr*tr^2/Ir-kr*b^2/Ip-kr/ms;...
        kf*t_f*tr/Ir+kf*a*b/Ip-kf/ms,-kf*t_f*tr/Ir+kf*a*b/Ip-kf/ms,kr*tr^2/Ir-kr*b^2/Ip-kr/ms,-kr*tr^2/Ir-kr*b^2/Ip-kr/ms];

D11 = [-t_f^2*Kv/Ir+a^2*Kv/Ip+Kv/ms,t_f^2*Kv/Ir+a^2*Kv/Ip+Kv/ms,-t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms,t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms;...
        t_f^2*Kv/Ir+a^2*Kv/Ip+Kv/ms,-t_f^2*Kv/Ir+a^2*Kv/Ip+Kv/ms,t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms,-t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms;...
       -t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms,t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms,-tr^2*Kv/Ir+b^2*Kv/Ip+Kv/ms,tr^2*Kv/Ir+b^2*Kv/Ip+Kv/ms;...
        t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms,-t_f*tr*Kv/Ir-a*b*Kv/Ip+Kv/ms,tr^2*Kv/Ir+b^2*Kv/Ip+Kv/ms,-tr^2*Kv/Ir+b^2*Kv/Ip+Kv/ms];

% Performance matrices for the 4 accelerations
C =[C11 zeros(4,11)];   % 
D =  zeros(4,4);    
F = zeros(4,4); 

% System's output
C12 = [eye(4);C11];
Cy  = [C12 zeros(8,11)];
Dyu = [zeros(4,4);D11];
Dyw = zeros(8,4);

% %%filter for the acceleration
% Fiso=1.0*tf([81.89 796.6 1937 0.1446],[1 80 2264 7172 21196]);% 
% Wa=append(Fiso,Fiso,Fiso,Fiso);
% [Ak Bk Ck Dk]=ssdata(Wa);

% Augmented matrices
Ag=A; %[A zeros(15,16);Bk*C  Ak];
Bgu=B;%[B; Bk*D];
Bgw=E;%[E;Bk*F];
Cgz=C;%[Dk*C  Ck];
Dgzu=D;% Dk*D; % 
Dgzw=F;%Dk*F;
% LQR 
X=sdpvar(4,4);% X is a scalar
Y=sdpvar(4,15);
P=sdpvar(15,15);
gamma=sdpvar(1);
G=[zeros(1,6) 1 zeros(1,4) -1 zeros(1,3);
    zeros(1,6) 1 zeros(1,5) -1 zeros(1,2);
    zeros(1,6) 1 zeros(1,6) -1 zeros(1,1);
    zeros(1,6) 1 zeros(1,7) -1];
Na=1/2;% 
G1=(Ag*P-Bgu*Y+P*Ag'-Y'*Bgu'<=0);
G2=([X  Y;...    
    Y'  P]>=0);
G3=(P>=0);
G4=((trace(-2*Kv*G'*Na*Y)+ Ra*trace(X))<=gamma);% trace(X)=X
G5=(gamma>=0);
 eta=10^5*sqrt(1);% 
Gri=([((Ag*P -Bgu*Y)'+(Ag*P-Bgu*Y))         Bgw             ((Cgz*P-Dgzu*Y))'; ...
     Bgw'                                 -eta*eye(4)           Dgzw'; ...
     (Cgz*P-Dgzu*Y)                        Dgzw           -eta*eye(4)]<=0);%
 ops=sdpsettings('solver','mosek');% 
 S=G1+G2+G3+G4+G5;
solvesdp(S, gamma ,ops);


 Kcl=double(Y)/double(P);
save('Kfull.mat', 'Kcl');
 tsim=150;
 sim('LQfullcar')

% 
