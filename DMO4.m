%������������ (��������������) ����� ���������� NOE

clear all; clc;close all;

start_time=clock;

width_x=5;
length_y=5;

m=1;j1=30000;L=1500; % measurement - samples - realizations

RazmX=2; % number of components of the state vector
 
X=zeros(RazmX,1); %�������������� �������� state vector
 
randn('state',0); rand('state',0); %��������� ������� ��������� ����� � ��������� ���������

b0=50; Px=diag([b0^2 b0^2]); % prior belief

P_NOE=zeros(RazmX,RazmX,m);P_NOE1=P_NOE;G_NOE1=P_NOE1; %P_NOE covariance matrix G_ effective
x1=zeros(L,1); x2=zeros(L,1); % coordinates for different realizations

xj1=zeros(1,j1); xj2=zeros(1,j1);% samples coordinates

% почему 7б а не 4?
Sij1=zeros(m,j1); Sij2=zeros(m,j1); Sij3=zeros(m,j1); Sij4=zeros(m,j1); Sij5=zeros(m,j1); Sij6=zeros(m,j1); 
Sij=zeros(4,j1);

% ������ ������� ������� � ����
SaveWorkSpace=0;

%���������� ���������� % coordinates of landmarks
Or1=[0,0]; %m
Or2=[width_x,length_y]; %m
Or3=[0,length_y]; %
Or4=[width_x,0]; %m

mX=750;h1=0; % что такой mX? = real position 

%���� ���������� что это? r1 RMS noise
r1=0.005*sqrt(((Or1(1)-mX)^2+(Or1(2)-mX)^2+h1^2)); r2=1*pi/180;
r3=r1; r4=r2; r5=r1;r6=r2;

sqrt(((Or1(1)-mX)^2+(Or1(2)-mX)^2+h1^2))

R_1=1.0*r1^2; 

R=diag(repmat(R_1,4,1));

R1=R^(-1);


Sd=sqrt(0.7*r1^2); Sp=sqrt(0.7*r2^2);

%����������� ���������
C1d=Sd*randn(1,j1);C1p=Sp*randn(1,j1);
Ci1d=Sd*randn(1,L);Ci1p=Sp*randn(1,L);

C2d=Sd*randn(1,j1);C2p=Sp*randn(1,j1);
Ci2d=Sd*randn(1,L);Ci2p=Sp*randn(1,L);

C3d=Sd*randn(1,j1);C3p=Sp*randn(1,j1);
Ci3d=Sd*randn(1,L);Ci3p=Sp*randn(1,L);

C4d=Sd*randn(1,j1);C4p=Sp*randn(1,j1);
Ci4d=Sd*randn(1,L);Ci4p=Sp*randn(1,L);


%������ ������ ����������
S_Or=0;
Er11=S_Or*randn(1,j1);Er12=S_Or*randn(1,j1);
Ed11=S_Or*randn(1,L);Ed12=S_Or*randn(1,L);

Er21=S_Or*randn(1,j1);Er22=S_Or*randn(1,j1);
Ed21=S_Or*randn(1,L);Ed22=S_Or*randn(1,L);

Er31=S_Or*randn(1,j1);Er32=S_Or*randn(1,j1);
Ed31=S_Or*randn(1,L);Ed32=S_Or*randn(1,L);

Er41=S_Or*randn(1,j1);Er42=S_Or*randn(1,j1);
Ed41=S_Or*randn(1,L);Ed42=S_Or*randn(1,L);

%�������� �������� ��������� �������
x1=mX+b0*randn(1,L); % gaussian distribution of realizations' first coordinate 
x2=mX+b0*randn(1,L); % gaussian distribution of realizations' second coordinate

x=[x1;x2];

% white noise - почему 8?
v1=sqrt(R(1,1))*randn(m,L);
v2=sqrt(R(2,2))*randn(m,L);
v3=sqrt(R(1,1))*randn(m,L);
v4=sqrt(R(2,2))*randn(m,L);
v5=sqrt(R(1,1))*randn(m,L);
v6=sqrt(R(2,2))*randn(m,L);
v7=sqrt(R(1,1))*randn(m,L);
v8=sqrt(R(2,2))*randn(m,L);
MaxErr=0; NumberErr=0;

mD=1/(10*60);mP=mD;

%vfnhbws lbyfvbrb b gjhj;lf.ob[ �����

dt=3*60; V0=0.15*0.514;

lc=1/(2*3600); d=exp(-lc*dt);

Fd=[1 0  dt    0 ;
    0 1  0    dt ;
    0 0  d     0 ;
    0 0  0     d     ];

mV=V0; d1=mV*sqrt(1-d^2); 

Gd=[0   0  ;
    0   0  ;
    d1  0  ;
    0  d1  ];

P2=[zeros(2,2)  zeros(2,2); zeros(2,2) diag([V0^2 V0^2])];GozX2=P2;
P3=zeros(2,2);G3=P3;

 %������������ ������� ������ �����-����� � ��������������� ������� ���������

        xj1=mX+b0*randn(1,j1); % samples distribution of the first coordinates 
        xj2=mX+b0*randn(1,j1); % samples distribution of the second coordinates
                     
        % Range_to_each landmark modeling for every samples
        for j=1:j1

                      
            Sij1(:,j)=sqrt( (xj1(1,j)- (Or1(1)+Er11(j)))^2 + (xj2(1,j)-(Or1(2)+Er12(j)) )^2 + h1^2 ) + C1d(j); 
            
            Sij3(:,j)=sqrt( (xj1(1,j)- (Or2(1)+Er21(j)))^2 + (xj2(1,j)-(Or2(2)+Er22(j)) )^2 + h1^2 ) + C2d(j); 
            
            Sij5(:,j)=sqrt( (xj1(1,j)- (Or3(1)+Er31(j)))^2 + (xj2(1,j)-(Or3(2)+Er32(j)) )^2 + h1^2 ) + C3d(j); 
            
            Sij7(:,j)=sqrt( (xj1(1,j)- (Or4(1)+Er41(j)))^2 + (xj2(1,j)-(Or4(2)+Er42(j)) )^2 + h1^2 ) + C4d(j); 
        
       
        end 
%
        Sij=[Sij1(1,:);Sij3(1,:);Sij5(1,:);Sij7(1,:)];
        
 %       R1=R1(1:4,1:4);
%         Sij=[Sij1(1,:);Sij2(1,:)];
        xj=[xj1;xj2];
       

for l=1:L %range for every realization


   y1=sqrt((x1(l)-(Or1(1)+Ed11(l)))^2+(x2(l)-(Or1(2)+Ed12(l)))^2 + h1^2) + Ci1d(l) + v1(:,l); %������������ ���������
   
   y3=sqrt((x1(l)-(Or2(1)+Ed21(l)))^2+(x2(l)-(Or2(2)+Ed22(l)))^2 + h1^2) + Ci2d(l) + v3(:,l); %������������ ���������
   
   y5=sqrt((x1(l)-(Or3(1)+Ed31(l)))^2+(x2(l)-(Or3(2)+Ed32(l)))^2 + h1^2) + Ci3d(l) + v5(:,l); %������������ ���������
   
   y7=sqrt((x1(l)-(Or4(1)+Ed41(l)))^2+(x2(l)-(Or4(2)+Ed42(l)))^2 + h1^2) + Ci4d(l) + v7(:,l); %������������ ���������
           
       
   
  y=[y1;y3;y5;y7];
 % y=[y3;y4;y5;y6];

  for k=1:m 
        
        Exp=0; Exp1=0; Exp2=0; Exp3=0;
        
        %����� �����-�����
       
        for h=1:j1
            
             Exp(h)=exp( -0.5*(y-Sij(:,h))'*R1*(y-Sij(:,h)) );
             
             Exp21=xj(:,h)*Exp(h);
                 
             Exp31=xj(:,h)*xj(:,h)'*Exp(h);
          
             Exp1=(1/j1)*Exp(h)+Exp1;
            
             Exp2=Exp21+Exp2;
                       
             Exp3=Exp31+Exp3;                    
        end

         OzX = (1/j1)*Exp2/(Exp1); %������
          
         P = (1/j1)*Exp3/Exp1 -OzX*OzX'; % ��������� ����������
                 
         GozX = (x(:,l)- OzX)*(x(:,l)- OzX)'; %�������������� ��������     
        
            if MaxErr < GozX
                MaxErr=GozX; NumberErr=l;
            end
        
            P2(1:2,1:2)=P; GozX2(1:2,1:2)=GozX;
            
        % P1=Fd*P2*Fd' + Gd*Gd';
         
        % GozX1=Fd*GozX2*Fd' + Gd*Gd';
            
            
 P_NOE1(:,:,k)=P2(1:2,1:2)+P_NOE1(:,:,k); %������������ ��� ������������ ����������
 
 P3=P+P3;
 
 G_NOE1(:,:,k)=GozX2(1:2,1:2)+G_NOE1(:,:,k);%������������ ��� ������������ ����������
 
 G3=GozX+G3;

      
end   
end

 P_NOE1=P_NOE1/L; %���������
 
 G_NOE1=G_NOE1/L;%��������������
 
 P3=P3/L;
 G3=G3/L;
 
dFi_OzG(1:m)=sqrt( G_NOE1(1,1,:) )
dLa_OzG(1:m)=sqrt( G_NOE1(2,2,:) )

dFi_OzP(1:m)=sqrt( P_NOE1(1,1,:) )
dLa_OzP(1:m)=sqrt( P_NOE1(2,2,:) )

OzX
x1(l)
x2(l)

% figure(3);plot(1:m,dFi_OzP,1:m,dFi_OzG,'LineWidth',2);grid %,'Color','red'
% figure(4);plot(1:m,dLa_OzP,1:m,dLa_OzG,'LineWidth',2);grid
%  
end_time=clock;

%start_time-end_time;

if SaveWorkSpace==1 save NOE1400TO.mat; end;
