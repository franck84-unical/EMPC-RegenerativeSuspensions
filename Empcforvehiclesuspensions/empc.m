clear all;

%% vehicle parameters
a1=1.4; % distanza tra la punta e il centro di massa [m]
a2=1.47; % distanza tra la coda e il centro di massa [m]
b1=0.7; % distanza tra la fiancata sinistra e il centro di massa [m]
b2=0.75; % distanza tra la fiancata destra e il centro di massa [m]
m=1274; %  Sprung mass[Kg] 
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

%% generator parameters

La=0.003;
Ka=48;
Ra=6;

%% Quarter car model parameters

ks=kf;  % rigidezza della sospensione
kt=ktf; % rigidezza del pneumatico
ms=m/4; % massa sospesa
mu=mf; % massa non sospesa

%% Road signal parameters
 n00=0.033;
n0=0.1;% reference spatial frequency, 
Gn0=80*10^-6;% this value is for a C road
pi=180;
alpha=-2*pi*n00;
beta=2*pi*n0*sqrt(Gn0);
u=70*1000/3600;%speed
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
Ts=10^-1;%Sampling time, the controller gain K has been designed offline considering this value
Ad=expm(A*Ts);
Bd=A^(-1)*(expm(A*Ts)-eye(5))*B;
Ed=A^(-1)*(expm(A*Ts)-eye(5))*E;% 
Cyd=Cy*A^(-1)*(expm(A*Ts)-eye(5));
G=[0 0 1 0 -1];%matrix expressing suspension's velocity w.r.t states
Na=1;%2*1/2
 n=5;m=1; 
 M=100;% the time steps
xreal=zeros(n,M);%real states

NT=15;   %prediction horizon
 
Verif=zeros(6*(NT+1),M);% I put inside the different solutions Z(t)
Pe1=zeros(1,M);% I put all the powers
w=zeros(1,M);% I put all the disturbances
sd1=zeros(1,M);% I put all the suspension deflections
td1=zeros(1,M);% I put all the tire deflections
acc1=zeros(1,M);% I put all the accelerations
u=zeros(1,M);%I put all the control inputs u=K*x
U=zeros(1,M);%I put all the control inputs U=v+K*x
EXIT=zeros(1,M);%I put all the exit values of fmincon

load('Wd.mat');% 100 values of the disturbance w variance 1
Sigma=1;%variance of the disturbance
%% the discrete time gain computed offline
load('K.mat');

%% objective function
H= [ G'*Na*K*Ka + K'*Ra*K            (G'*Na*Ka +K'*Ra);
     Ra*K                             Ra   ];
H=(H+H')/2; % just to make it symmetric, it doesn't impact much the results. 
HT=H;
for i=1:NT-1
HT=blkdiag(HT,H);%global objective function
end
phi=eye(5);
Hplus=[phi       zeros(5,1);
       zeros(1,5)    0];%recursive feasibility condition
HT=blkdiag(HT,Hplus);
 Xtilde=zeros(5,5);
F=trace((G'*Na*Ka*K+K'*Ra*K)*Xtilde);
for i=1:NT-1
Xtilde=(Ad+Bd*K)*Xtilde*(Ad+Bd*K)'+Ed*Sigma*Ed';
F=F+trace((G'*Na*Ka*K+K'*Ra*K)*Xtilde);
end
F=F+trace(phi*Xtilde);
%% equality constraints  Aeq*Z=Beq
aeq=[ Ad+Bd*K   Bd  -eye(5,5)];%5*11
Aeq=zeros( 5*(NT+1)+5, 6*(NT+1));%
Aeq(1:5,:)=[eye(5), zeros(5,6*(NT+1)-5)];
 for k=0:(NT-1)
        Aeq((6+5*k):(10+5*k),:)=[zeros(5,6*k),aeq,zeros(5,(6*(NT+1)-6*k-11))];
 end
  Aeq(5*(NT+1)+1:5*(NT+1)+5,:)=[ zeros(5,6*(NT+1)-11), aeq];
 Beq=zeros(5*(NT+1)+5,1);
Beq(1:5)=xreal(:,1);%xequilibrium at time t=1, initial

  T=zeros(1,M);%time vector

 %% initial guess 
Z0=0*ones(6*(NT+1),1);% initial guess
  

 %% MPC loop
    for t=1:M

%% inequalities constraints
% constraint on the input u
umax=25;
aunu=[K 1];
gx4=umax;

 % constant constraint tightening on the acceleration
    Fxacc=[0 1 0 -1 0];
 aunacc=(1/ms)*(-ks* [Fxacc  0]  + Ka*[zeros(1,5) 1] );
 aun1=aunacc;
 %constant constraint on tire deflection and suspension deflection
   Fx2=[-1 0 0 1 0]; 
  aun2=[Fx2      0]; % for tire defletion
   Fx3=[0 1 0 -1 0]; aun3=[Fx3  0]; % for suspension deflection
  % time varying constraint tightening on tire deflection
   delta=0; deltasd=0; X=0;
gx1=3; 
gx2=0.07-delta;%2nd approach
 gx3=0.07-deltasd;
%% to apply the constraints, just comment and uncomment the constraints you want. 
%AUN=[aun1;-aun1;aun2;-aun2;aun3;-aun3]; BUN=[gx1;gx1;gx2;gx2;gx3;gx3];%1st approach
%AUN=[aun1;-aun1;aunu;-aunu]; BUN=[gx1;gx1;gx4;gx4];%4th approach
AUN=[aun2;-aun2;aun3;-aun3]; BUN=[gx2;gx2;gx3;gx3];%2nd approach
%AUN=[aun1;-aun1]; BUN=[gx1;gx1];%3rd approach
for i=1:NT


%aun=[aun1;-aun1;aun2;-aun2;aun3;-aun3];%1st approach
aun=[aun2;-aun2;aun3;-aun3];%2nd approach
%aun=[aun1;-aun1];%3rd approach
%aun=[aun1;-aun1;aunu;-aunu]; %4th approach

%bun=[gx1;gx1;gx2;gx2;gx3;gx3];%1st approach
bun=[gx2;gx2;gx3;gx3];%2nd approach
%bun=[gx1;gx1];%3rd approach
%bun=[gx1;gx1;gx4;gx4];%4th approach
AUN=blkdiag(AUN,aun);
 X=(Ad+Bd*K)*X*(Ad+Bd*K)'+Ed*Sigma*Ed';
  delta=sqrt(trace(Fx2'*Fx2*X));%update of delta,tire deflection,2nd approach
  deltasd=sqrt(trace(Fx3'*Fx3*X));%update of deltasd, suspension deflection,2nd approach
  gx2=0.07-delta;%2nd approach
  gx3=0.07-deltasd;

BUN=[BUN;bun];
end  

M=[zeros(5,5)      G'*Na*Ka;
     zeros(1,5)      Ra      ];% M is 6*6 square non symmetric
 M=(M+M')/2;%we symmetrize


%% final objective function
fun = @(Z) ( Z'*HT*Z +F);%removing F doesn't change the result. 
  NONLCON=@nonlcon;%this contains the terminal constraint
     %AUN=[];BUN=[];% 
      Lb=[];Ub=[];% 
   MyValue = 1e+10;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',MyValue);
           [Z,FVAL,EXITFLAG] = fmincon(fun,Z0,AUN,BUN,Aeq,Beq,Lb,Ub,NONLCON,options);%
              EXITFLAG% a digit which indicates the feasibility of the problem
    EXIT(1,t)=EXITFLAG;% I store the digits in the vector EXIT

       %% real x(t) and electrical power
  U(:,t)=Z(6)+K*xreal(:,t);% I store the control inputs
  Pe1(1,t)= xreal(:,t)'*G'*Na*Ka*U(:,t)+U(:,t)'*Ra*U(:,t);% compute the electrical power
  disp(t);
  Verif(:,t)=Z; %I store the solution vector Z(t) 
  T(1,t)=t;
  xreal(:,t+1)=(Ad+Bd*K)*xreal(:,t)+Bd*Z(6)+ Ed*w(1,t);%next states 
  u(:,t)=K*xreal(:,t);%I store the control inputs
  disp(Z(6)+K*xreal(:,t));
  disp(Pe1(1,t));
    sd1(1,t)=xreal(2,t)-xreal(4,t);  td1(1,t)=xreal(4,t)-xreal(1,t);
  acc1(1,t)=(1/ms)*(-ks*(xreal(2,t)-xreal(4,t))+Ka*U(:,t));
  if (EXITFLAG<0)
       break;
   end
  %% the new initial guess and new initial states
  Z0=Z;
  Beq(1:5)=xreal(:,t+1);%xequilibrium at time t
%% method to remove unfeasibility

    end

    mean(Pe1)
    %% comparison plots
load('acclq.mat');% I load the LQ results
load('sdlq.mat');
load('tdlq.mat');
load('Pelq.mat');
load('ulq.mat');
subplot(5,1,1)%power
hold on; grid on;
plot(T,Pe,'b') ;
plot(T,Pe1,'r');
ylabel('Power[W]','Interpreter','latex','FontSize',12);

subplot(5,1,2)%current
hold on; grid on;
plot(T,ulq,'b') ;
plot(T,U,'r');
ylabel('u[A]','Interpreter','latex','FontSize',12);

subplot(5,1,3)%acceleration
hold on; grid on;
plot(T,acc,'b') ;
plot(T,acc1,'r');
set(gca,'FontSize',10,'FontWeight','Bold')
ylabel('$\ddot{z}_{s} [m/s^2]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
%ylabel('$\ddot{z}_{1} [m/s^2]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
subplot(5,1,4)%tire deflection
hold on; grid on;
plot(T,td,'b') ;
plot(T,td1,'r');
set(gca,'FontSize',10,'FontWeight','Bold')
ylabel('${z}_{u} -{z}_{r} [m]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
%ylabel('${z}_{2} -{z}_{3} [m]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
subplot(5,1,5)%suspension deflection 
hold on; grid on;
plot(T,sd,'b') ;
plot(T,sd1,'r');
legend('LQR','Constrained EMPC 2','FontSize',10,'FontWeight','Bold')
%legend('LQR','Linear EMPC','FontSize',10,'FontWeight','Bold')
set(gca,'FontSize',10,'FontWeight','Bold')
xlabel('time[s]', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
ylabel('${z}_{s} -{z}_{u} [m]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
%ylabel('${z}_{1} -{z}_{2} [m]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');