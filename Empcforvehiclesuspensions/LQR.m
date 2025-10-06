%% LQ control in discrete time 
clear all;
%% vehicle parameters
a1=1.4; % distanza tra la punta e il centro di massa [m]
a2=1.47; % distanza tra la coda e il centro di massa [m]
b1=0.7; % distanza tra la fiancata sinistra e il centro di massa [m]
b2=0.75; % distanza tra la fiancata destra e il centro di massa [m]
m=1274; %  Sprung mass  [Kg] 
mf=35.5; % 35.5 front wheel mass [Kg]
mr=35.5; % 35.5 rear wheel mass [Kg]
cf=1140; % smorzamento degli ammortizzatori anteriori [N*s/m]
cr=1250; % smorzamento degli ammortizzatori posteriori [N*s/m]
kf=27000; % 27000 Front suspension stiffness[N/m]
kr=30000; % 30000 Rear suspension stiffness [N/m]
Ix=1523; % 1523 pitch inertia momento di inerzia attorno all'asse di becheggio [Kg*m^2]
Iy=606.1; % 606.1 roll inertia [Kg*m^2]
kR=25000; % rigidezza torsionale della barra antirollio [N*m/rad]
ktf=228000; % 228000 Front tyre stiffness [N/m]
ktr=228000; % 228000 Rear tyre stiffness [N/m]
tsim=150;

% generator parameters

La=0.003;
Ka=48;
Ra=6;

% Quarter car model parameters

ks=kf;  % rigidezza della sospensione
kt=ktf; % rigidezza del pneumatico
ms=m/4; % massa sospesa
mu=mf; % massa non sospesa

%% Road signal parameters
 n00=0.033;
n0=0.1;% reference spatial frequency page 5, 
Gn0=80*10^-6;%  this value is for a C road
pi=180;
alpha=-2*pi*n00;
beta=2*pi*n0*sqrt(Gn0);
u=70*1000/3600;
%% Continuous time model
 A=[alpha*u zeros(1,4);
    zeros(1,2) 1 zeros(1,2);
    0 -ks/ms 0 ks/ms 0;
    zeros(1,4)  1 ;
    kt/mu ks/mu 0 -(ks+kt)/mu 0];
 
B=[0
   0
   Ka/ms
   0
   -Ka/mu];

E=[beta*sqrt(u)
    0
   0
   0
   0];

Cy=[0 0 1 0 -1];
Dyu=0;
Fyw=0;


C=[0 -ks/ms 0  ks/ms 0 ];%Only the acceleration

Dzu=Ka/ms;

Fzw=0;

 %% discrete time model
Ts=10^-1;%sampling time
Ad=expm(A*Ts);
Bd=A^(-1)*(expm(A*Ts)-eye(5))*B;
Ed=A^(-1)*(expm(A*Ts)-eye(5))*E;% 
Cyd=Cy*A^(-1)*(expm(A*Ts)-eye(5));
G=[0 0 1 0 -1];
Na=1;
 M=[zeros(5,5)      G'*Na*Ka;
     zeros(1,5)      Ra      ];% M is 6*6 square non symmetric
 M=(M+M')/2;%M is 6*6 symmetric and positive semi definite.  
 M11=zeros(5);
 M12=[0;0;24;0;-24];
 M22=6;

Sigma=1;% variance of the signal
load('Wd.mat');% 100 values of the signal w of variance 1 

M=10^2;%time steps

Pe=zeros(1,M);%
 sd=zeros(1,M);% I put all the suspension deflection
td=zeros(1,M);% I put all the tire deflections
acc=zeros(1,M);% I put all the accelerations
T=zeros(1,M);
i=1;

 eta=10^3;   
%% LMI
PI=sdpvar(5,5);
Y=sdpvar(1,5);
Z=sdpvar(1,1);
gamma=sdpvar(1,1);
N1=sdpvar(1,1);   X=sdpvar(1,1);
G1=([PI                       PI*Ad'+Y'*Bd';...    
    (PI*Ad'+Y'*Bd')'        PI-Ed*Sigma*Ed']>=0);%this is an  equality constraint originally.
I1=M22^(1/2)*Y;
G2=([Z  I1;...    
    I1'  PI]>=0);
G3=(trace(M11*PI+Y'*M12'+M12*Y)+Z<=gamma);
G4=(gamma>=-eta);%-10^3 
G5=(PI>=0);

 S=G1+G2+G3+G4+G5;
ops=sdpsettings('solver','mosek');% 

solvesdp(S,gamma,ops)
 K=double(Y)/double(PI);
 if (abs(eig(Ad+Bd*K))<1) % stability condition
    disp(K);
    %gamma
    vecpi=inv(eye(25)-kron((Ad+Bd*K),(Ad+Bd*K)))*reshape(Ed*Sigma*Ed',25,1);
    PI=reshape(vecpi,5,5);
    Power=trace((M12*K+K'*M12'+K'*M22*K)*PI)%maximum possible mean power 
 end
% 
 %% computation of the real power, control input, suspension, tire and acceleration at each time step
 n=5;m=1; 
 xreal=zeros(n,M);
 xreal(:,1)=0*ones(5,1); %initial states
 ulq=zeros(1,M);%initial control input
 for t=1:M
  Pe(1,t)=(Ka*xreal(:,t)'*G'*Na*(0+K*xreal(:,t))+(0+K*xreal(:,t))'*Ra*(0+K*xreal(:,t)));
 % disp(t);
  T(1,t)=t;
  xreal(:,t+1)=(Ad+Bd*K)*xreal(:,t)+Bd*0+Ed*w(1,t);% next real states
  ulq(:,t)=K*xreal(:,t);% control input update
  sd(1,t)=xreal(2,t)-xreal(4,t);  td(1,t)=xreal(4,t)-xreal(1,t);
  acc(1,t)=(1/ms)*(-ks*(xreal(2,t)-xreal(4,t))+Ka*K*xreal(:,t));
 end

mean(Pe(1,:))
%% I save the control gain and the obtained results
save('K.mat', 'K');
save('Pelq.mat', 'Pe');
save('ulq.mat', 'ulq');
save('acclq.mat', 'acc');
save('tdlq.mat', 'td');
save('sdlq.mat', 'sd');