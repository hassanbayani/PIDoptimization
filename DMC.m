% This program performs DMC controller for MIMO systems

function [ys,y,ym,u,du]=DMC(plant,delayp,NU,NY,g,par,alfa,gama,q,D,T)
% N may be convert to a (NU X NY) matrix differ for every input & output

M=par(1); n1=par(2); n2=par(3); N=par(4); Ts=par(5); sw=par(6);
P=n2-n1+1;

plant=D*plant;
set(plant,'iodelay',delayp);

plantd=c2d(plant,Ts,'zoh');
delp=plantd.iodelay;

for i=1:NY;
    for j=1:NU,
        B=plantd.num{i,j};
        A=plantd.den{i,j};
        Bp(i,j,1:length(B)+delp(i,j))=[zeros(1,delp(i,j)) B];
        Ap(i,j,1:length(B)+delp(i,j))=[A zeros(1,delp(i,j))];
        nbp(i,j)=length(Bp);
        nap(i,j)=length(Ap);
    end,
end,

switch sw,
    case 1,
        for i=1:NY,
            ys(i,:)=[zeros(1,(i-1)*fix(0.5*T/NY)) ones(1,fix(0.5*T)) zeros(1,fix(0.5*T)-(i-1)*fix(0.5*T/NY))];
        end,
    case 2,
        for i=1:NY,
            ys(i,:)=sin(2*pi*0.002*(0:Ts:Ts*(T-1))/Ts+(i-1)*pi/NY);
        end,
    case 3,
        for i=1:NY,
            ys(i,:)=square(2*pi*0.002*(0:Ts:Ts*(T-1))/Ts+(i-1)*pi/NY);
        end,
end,

% Constructs dynamic matrix of the system

G=zeros(NY*P,NU*M);

for i=1:NY
    for j=1:NU
        A1=zeros(M,1);
        yy=g(:,i,j);
        A1(1)=yy(n1);
        A2=yy(n1:n2);
        G((i-1)*P+1:i*P,(j-1)*M+1:j*M)=Toeplitz(A2,A1);
    end
end

GN=zeros(NY*P,(N-1)*NU); % Hankel  P*NY X (N-1)*NU
for i=1:NY
    for j=1:NU
        R=zeros(N-1,1);
        C=zeros(n2,1);
        yy=g(:,i,j);
        C=yy(3:n2+2);
        R(1:N-n2)=yy(n2+2:N+1);
        GN((i-1)*n2+1:i*n2,(j-1)*(N-1)+1:j*(N-1))=Hankel(C,R);
    end,
end,

% To calculate Ypast

for i=1:NY,
    for j=1:NU,
        gdc(i,j)=g(N+2,i,j);
        g1(i,j)=g(2,i,j);
    end,
end,

gN=zeros(NY*n2,NU*n2);

for i=1:NY,
    for j=1:NU,
        for l=1:n2,
            gN((i-1)*n2+l,(j-1)*n2+l)=gdc(i,j);
        end,
    end,
end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate KDMC and Define Control Matrix (Q,R)

qq=[];
for i=1:NY,
    qq=[qq;q(i,i)*ones(n2-n1+1,1)];
end,
Q=diag(qq);

rr=[];
for j=1:NU,
    rr=[rr;gdc(j,j)*gama(j,j)*ones(M,1)];
end,
R=diag(rr);

KDMC=inv(G'*Q*G+R)*G'*Q;

% size(KDMC) NU.M*NY.P
u(1:NU,1)=zeros(NU,1);
ut(1:NU,1)=zeros(NU,1);
uN=zeros(NU*n2,1);
uNi=zeros(NU*N,1);
duN=zeros(NU*(N-1),1);
y=zeros(NY,1); % Yplant
ypast=zeros(NY*n2,1);
ym(1:NY,1)=zeros(NY,1);
dd(1:NY,1)=zeros(NY,1);

% Fixed Control Matrix for Model

for k=1:T,

    % Y desired Calculation for each step

    for i=1:NY,
        ydd(i,k)=y(i,k);
        for l=1:n2,
            ydd(i,k+l)=alfa(i,i)*ydd(i,k+l-1)+(1-alfa(i,i))*ys(i,k);
        end,
        yd((i-1)*P+1:i*P,1)=ydd(i,k+n1:k+n2)';
    end,

    %Ypast Calculation for model********************
    ypast=GN*duN+gN*uN;
    
    for i=1:NY,
        yp((i-1)*P+1:i*P,1)=ypast((i-1)*n2+n1:(i-1)*n2+n2)';
    end,

    for i=1:NY,
        d((i-1)*P+1:i*P,1)=dd(i,k)*ones(P,1);
    end,

    duu =KDMC*(yd-d-yp);

    for j=1:NU,
        du(j,k)=duu((j-1)*M+1);
    end,

    %updating u
    for j=1:NU,
        u(j,k)=ut(j)+du(j,k);
        ut(j)=u(j,k);
    end,
    
    for j=1:NU,
        for l=1:N-1,
            uNi((j-1)*N+l)=uNi((j-1)*N+l+1);
        end,
    end,

    for j=1:NU,
        uNi((j-1)*N+N)=u(j,k);
    end,
    for j=1:NU,
        for l=1:n2,
            uN((j-1)*n2+l)=uNi((j-1)*N+l);
        end,
    end,
    
    for j=1:NU,
        for l=N-1:-1:2,
            duN((j-1)*(N-1)+l)=duN((j-1)*(N-1)+l-1);
        end,
    end,

    for j=1:NU,
        duN((j-1)*(N-1)+1)=du(j,k);
    end,

    for i=1:NY,
        for j=1:NU,
            yt(i,j,k+1)=0;
            for n=2:nap(i,j),
                if k+2-n > 0,
                    yt(i,j,k+1)=yt(i,j,k+1)-Ap(i,j,n)*yt(i,j,k+2-n);
                end,
            end,
            for n=2:nbp(i,j),
                if k+2-n > 0,
                    yt(i,j,k+1)=yt(i,j,k+1)+Bp(i,j,n)*u(j,k+2-n);
                end,
            end,
        end,
        y(i,k+1)=0;
        for j=1:NU,
            y(i,k+1)=y(i,k+1)+yt(i,j,k+1);
        end,
    end,

    for i=1:NY,
        ym(i,k+1)=ypast((i-1)*n2+1);
    end,
    
    ym(:,k+1)=ym(:,k+1)+g1*du(:,k);

    for i=1:NY,
        dd(i,k+1)=y(i,k+1)-ym(i,k+1);
    end,
end,

for j=1:NU,
    u(j,T-1)=u(j,T-2);
    du(j,T-1)=du(j,T-2);
end,

NN=max(NY,NU);
figure(1)
for i=1:NY,
    subplot(4,NN,i), plot((1:T)*Ts,y(i,1:T),'k--',(1:T)*Ts,ys(i,1:T),'k-'), ylabel('ys, y'), grid
end,
for i=1:NY,
    subplot(4,NN,NN+i), plot((1:T)*Ts,y(i,1:T),'k--'), ylabel('ym'), grid
end,
for j=1:NU,
    subplot(4,NN,2*NN+j), plot((1:T)*Ts,u(j,1:T),'k--'), ylabel('u'), grid
end,
for j=1:NU,
    subplot(4,NN,3*NN+j), plot((1:T)*Ts,du(j,1:T),'k--'), ylabel('du'), grid
end,
%%%%%%%%%khodam neveshtam

% 
% for j=1:NU,
%     figure(j+1)
%     plot((1:T)*Ts,u(j,1:T),'k--'), ylabel('u'), grid
% end,