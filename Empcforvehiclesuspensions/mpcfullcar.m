clear all
clc;


yalmip('clear')

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

% generator parameters

La=0.003;
Kv=48;
Ka=48;
Ra=6;


c=100; % 
%% Road signal parameters
 n00=0.033;
n0=0.1;% reference spatial frequency page 5, 
Gn0=80*10^-6;%  this value is for a C road
pi=180;
alpha=-2*pi*n00;
beta=2*pi*n0*sqrt(Gn0);
u=70*1000/3600;

%% Full car model parameters
set_p=[ms;muf;mur;Ip;Ir;kf;kr;t_f;tr;ktf;ktr;bf;br;a;b];

set_a=[La;Kv;Ra];

% Submatrices
A12 = [t_f a 1; -t_f a 1; tr -b 1; -tr -b 1];
A21 = [-kf*t_f/Ir kf*t_f/Ir -kr*tr/Ir -kr*tr/Ir; -kf*a/Ip -kf*a/Ip kr*b/Ip kr*b/Ip; -kf/ms -kf/ms -kr/ms -kr/ms];
A31 = diag([kf/muf,kf/muf,kr/mur,kr/mur]);
A33 = -diag([ktf/muf,ktf/muf,ktr/mur,ktr/mur]);
% A31 = [kf/muf 1 1 1; 1 kf/muf 1 1; 1 1  kr/mur 1; 1 1 1 kr/mur];
% A33 = -[ktf/muf 1 1 1; 1 ktf/muf 1 1; 1 1 ktr/mur 1; 1 1 1 ktr/mur];
B21 = [t_f*Kv/Ir -t_f*Kv/Ir tr*Kv/Ir -tr*Kv/Ir; a*Kv/Ip a*Kv/Ip -b*Kv/Ip -b*Kv/Ip; Kv/ms Kv/ms Kv/ms Kv/ms];
B41 = -diag([Kv/muf,Kv/muf,Kv/mur,Kv/mur]);

% State Matrices
A = [zeros(4,4) A12 zeros(4,4) -eye(4); A21 zeros(3,3) zeros(3,4) zeros(3,4); zeros(4,4) zeros(4,3) zeros(4,4) eye(4); A31 zeros(4,3) A33 zeros(4,4)];
%A = [10^2*ones(4,4) A12 10*ones(4,4) -10*ones(4); A21 10*ones(3,3) ones(3,4) ones(3,4); ones(4,4) ones(4,3) ones(4,4) ones(4); A31 ones(4,3) A33 ones(4,4)];
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
G=[zeros(1,6) 1 zeros(1,4) -1 zeros(1,3);
    zeros(1,6) 1 zeros(1,5) -1 zeros(1,2);
    zeros(1,6) 1 zeros(1,6) -1 zeros(1,1);
    zeros(1,6) 1 zeros(1,7) -1];
 %% discrete time model
Ts=10^-1;%don't change, the K has been designed considering this value
sysc=ss(A,B,Cy,Dyu);
sysd=c2d(sysc,Ts);
sysc2=ss(A,E,Cy,Dyw);
sysd2=c2d(sysc2,Ts);
Ad=sysd.A;
Bd=sysd.B;
Cyd=sysd.C;
Ed=sysd2.B;
Na=1;
 n=15;m=1; 
 M=100;% time steps, 8
xreal=zeros(n,M);
%Prediction horizon
NT=15;   
 
Verif=zeros(19*(NT+1),M);% I put inside the different solutions Z(t)
Petot=zeros(1,M);% I put all the powers
Pe1=zeros(1,M);
Pe2=zeros(1,M);
Pe3=zeros(1,M);
Pe4=zeros(1,M);
w=zeros(4,M);% I put all the disturbances
sd1=zeros(1,M);% I put all the suspension deflection
td1=zeros(1,M);% I put all the tire deflection
sd2=zeros(1,M);% I put all the suspension deflection
td2=zeros(1,M);% I put all the tire deflection
acc1=zeros(1,M);% I put all the accelerations
acc2=zeros(1,M);% I put all the accelerations
 U=zeros(4,M);
  u=zeros(4,M);
EXIT=zeros(1,M);

Sigma=1;

load('Kfull');
K=Kcl;
%% objective function
H= [ G'*Na*Ka*K + K'*Ra*K            (G'*Na*Ka +K'*Ra);
     Ra*K                             Ra*eye(4)   ];
H=(H+H')/2; % doesn't impact the results
HT=H;
for i=1:NT-1
HT=blkdiag(HT,H);
end
%recursive feasibility
phi=eye(15);%% recursive feasibility
Hplus=[phi       zeros(15,4);
       zeros(4,15)    zeros(4,4)];
HT=blkdiag(HT,Hplus);
 Xtilde=zeros(15,15);
F=trace((G'*Na*Ka*K+K'*Ra*K)*Xtilde);
for i=1:NT-1
Xtilde=(Ad+Bd*K)*Xtilde*(Ad+Bd*K)'+Ed*Sigma*Ed';
F=F+trace((G'*Na*Ka*K+K'*Ra*K)*Xtilde);
end
F=F+trace(phi*Xtilde);
%% equality constraints  Aeq*Z=Beq
aeq=[ Ad+Bd*K   Bd  -eye(15)];%15*34
Aeq=zeros( 15*(NT+1)+15, 19*(NT+1));%
Aeq(1:15,:)=[eye(15), zeros(15,19*(NT+1)-15)];
 for k=0:(NT-1)%check this again
        Aeq((16+15*k):(30+15*k),:)=[zeros(15,19*k),aeq,zeros(15,(19*(NT+1)-19*k-34))];
 end
  Aeq(15*(NT+1)+1:15*(NT+1)+15,:)=[ zeros(15,19*(NT+1)-34), aeq];
 Beq=zeros(15*(NT+1)+15,1);
Beq(1:15)=xreal(:,1);%xequilibrium at time t=1, initial

  T=zeros(1,M);


 %% initial guess and initial states
Z0=0*ones(19*(NT+1),1);% initial guess
   xreal(:,1)=0*ones(15,1);   
T=zeros(1,M);
%  W=zeros(4,M);
load('Wd.mat');
w=0.3*w;
  W=[w;w;w;w];
 %% MPC loop
    for t=1:M

%% inequalities constraints
% constraint on the input u
umax=500;
aunu=[K eye(4)];
gx4=[umax;umax;umax;umax];
 
 % constant constraint tightening on the acceleration
%  acc(1,t)=C(2,:)*xreal(:,t) + D(2,:)*ulq(:,t) + F(2,:)*w(:,t);
  aunacc=[C(2,:) zeros(1,4)]; 
  aun1=aunacc;
%acc1(1,t)=(1/ms)*(-kf*xreal(1,t) -kf*xreal(2,t) -kr*xreal(3,t)-kr*xreal(4,t)+Ka*(U(1,t)+U(2,t)+U(3,t)+U(4,t))  );
 %constant constraint tightening on tire deflection
   aun2=[zeros(1,8) 1 zeros(1,10)]; 
 aun3=[0 1 zeros(1,17)]; % for suspension deflection
   gx2=0.1;%tire deflection
 gx3=0.1;%suspension deflection
gx1=5; %acceleration

% AUN=[aun2;-aun2;aun3;-aun3;aunu;-aunu]; BUN=[gx2;gx2;gx3;gx3;gx4;gx4];
AUN=[aun1;-aun1;aun2;-aun2;aun3;-aun3]; BUN=[gx1;gx1;gx2;gx2;gx3;gx3];
for i=1:NT
%aun=[aun2;-aun2;aun3;-aun3;aunu;-aunu];
aun=[aun1;-aun1;aun2;-aun2;aun3;-aun3];
AUN=blkdiag(AUN,aun);
gx2=0.1;
 gx3=0.1;
 gx1=5;
bun=[gx1;gx1;gx2;gx2;gx3;gx3];
BUN=[BUN;bun];
end  


%% big objective function
fun = @(Z)( Z'*HT*Z +F);%removing F doesn't change the result. 
  NONLCON=[];%@nonlconfullcar; %leads to infeasibilities
    AUN=[];BUN=[];% 
      Lb=[];Ub=[];% 
         %% initial conditions   
 MyValue = 1e+10;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',MyValue);
           [Z,FVAL,EXITFLAG] = fmincon(fun,Z0,AUN,BUN,Aeq,Beq,Lb,Ub,NONLCON,options);%
              EXITFLAG% a digit which indicates the feasibility of the problem
    EXIT(1,t)=EXITFLAG;% I store the digits in a vector

% if (t==63)
%     Z=Verif(:,62);
% end
       %% real x(t) and electrical power
         U(:,t)=Z(16:19)+K*xreal(:,t);
      % Petot(:,t)= xreal(:,t)'*G'*Na*Ka*U(:,t)+U(:,t)'*Ra*U(:,t);
      w=W(:,t);
 H1=[zeros(1,4) t_f a 1 zeros(1,8)]; U1=[zeros(1,11) 1 zeros(1,3)];
  H2=[zeros(1,4) -t_f a 1 zeros(1,8)]; U2=[zeros(1,12) 1 zeros(1,2)];
  H3=[zeros(1,4) tr -b 1 zeros(1,8)]; U3=[zeros(1,13) 1 zeros(1,1)];
  H4=[zeros(1,4) -tr -b 1 zeros(1,8)]; U4=[zeros(1,14) 1];
  A1=(1/Ir)*[-kf*t_f kf*t_f -kr*tr kr*tr zeros(1,11)];
  B1=(1/Ir)*[t_f*Ka -t_f*Ka tr*Ka -tr*Ka];%expression of dotdotphis w.r.t x and u
  A2=(1/Ip)*[-kf*a  -kf*a kr*b  kr*b  zeros(1,11)];
  B2=(1/Ip)*[a*Ka a*Ka -b*Ka -b*Ka];%expression of dotdotthetas w.r.t x and u
  %I compute the power at each wheel
   Pe1(1,t)=Kv*xreal(:,t)'*(H1-U1)'*Na*U(1,t) +  U(1,t)'*Ra*U(1,t);
        Pe2(1,t)=Kv*xreal(:,t)'*(H2-U2)'*Na*U(2,t) +  U(2,t)'*Ra*U(2,t);
           Pe3(1,t)=Kv*xreal(:,t)'*(H3-U3)'*Na*U(3,t) +  U(3,t)'*Ra*U(3,t);
              Pe4(1,t)=Kv*xreal(:,t)'*(H4-U4)'*Na*U(4,t) +  U(4,t)'*Ra*U(4,t);
   Petot(1,t)=Pe1(1,t)+Pe2(1,t)+Pe3(1,t)+Pe4(1,t);%total power
  disp(t);
  Verif(:,t)=Z; %I store Z(t) 
  T(1,t)=t;
     xreal(:,t+1)=(Ad+Bd*K)*xreal(:,t)+Bd*Z(16:19)+ Ed*w;%next states 
  u(:,t)=K*xreal(:,t);
    sd1(1,t)=Z(1);  td1(1,t)=Z(8);
    sd2(1,t)=Z(2);  td2(1,t)=Z(9);
  acc1(1,t)=(1/ms)*(-kf*xreal(1,t) -kf*xreal(2,t) -kr*xreal(3,t)-kr*xreal(4,t)+Ka*(U(1,t)+U(2,t)+U(3,t)+U(4,t))  );%correct this
  %the acceleration of the vehicle
  acc2(1,t)=-t_f*(A1*xreal(:,t)+B1*U(:,t)) + a*(A2*xreal(:,t)+B2*U(:,t)) + acc1(1,t);
  %expression of dotdotzs w.r.t x and u
  if (EXITFLAG<0)
       break;
  end
  %% the new initial guess and new initial states
  Z0=Z;
  Beq(1:15)=xreal(:,t+1);%xequilibrium at time t


    end
    mean(Petot)
% load('acclq.mat');
% load('sdlq.mat');
% load('tdlq.mat');
% load('Pelq.mat');
% load('ulq.mat');
% subplot(5,1,1)%power
% hold on; grid on;
% plot(T,Pe,'b') ;
% plot(T,Pe1,'r');
% ylabel('Power[W]','Interpreter','latex','FontSize',12);
% 
% subplot(5,1,2)%current
% hold on; grid on;
% plot(T,ulq,'b') ;
% plot(T,U,'r');
% ylabel('u[A]','Interpreter','latex','FontSize',12);
% 
% subplot(5,1,3)%acceleration
% hold on; grid on;
% plot(T,acc,'b') ;
% plot(T,acc1,'r');
% set(gca,'FontSize',10,'FontWeight','Bold')
% ylabel('$\ddot{z}_{s} [m/s^2]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
% 
% subplot(5,1,4)%tire deflection
% hold on; grid on;
% plot(T,td,'b') ;
% plot(T,td1,'r');
% set(gca,'FontSize',10,'FontWeight','Bold')
% ylabel('${z}_{u} -{z}_{r} [m]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
% 
% subplot(5,1,5)%suspension deflection 
% hold on; grid on;
% plot(T,sd,'b') ;
% plot(T,sd1,'r');
% legend('LQR','Constrained EMPC','FontSize',10,'FontWeight','Bold')
% set(gca,'FontSize',10,'FontWeight','Bold')
% xlabel('time[s]', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
% ylabel('${z}_{s} -{z}_{u} [m]$', 'Interpreter','latex','FontSize',12,'FontWeight','Bold');
