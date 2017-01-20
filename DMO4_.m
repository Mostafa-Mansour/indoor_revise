%������������ (��������������) ����� ���������� NOE

clear all; clc;close all;

start_time=clock;

width_x=5;
length_y=5;

m=1;j1=30000;L=100; % measurement - samples - realizations

RazmX=2; % number of components of the state vector
 
X=zeros(RazmX,m,L); %�������������� �������� state vector
 
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
r1=0.005*sqrt(((Or1(1)-mX)^2+(Or1(2)-mX)^2+h1^2)); 


R_1=1.0*r1^2; 

R=diag(repmat(0.25,4,1));

R1=R^(-1);


% Sd=sqrt(0.7*r1^2); 
Sd=0;

% modeling systematic for samples and realizations
C1d=Sd*rand(1,j1); %samples
Ci1d=Sd*rand(1,L); %realizations

C2d=Sd*rand(1,j1);
Ci2d=Sd*rand(1,L);

C3d=Sd*rand(1,j1);
Ci3d=Sd*rand(1,L);

C4d=Sd*rand(1,j1);
Ci4d=Sd*rand(1,L);




%�������� �������� ��������� �������
x1=0.25+(width_x-0.25)*rand(1,L); % gaussian distribution of realizations' first coordinate 
x2=0.25+(length_y-0.25)*rand(1,L); % gaussian distribution of realizations' second coordinate

x=[x1;x2];

% white noise - почему 8?
v1=sqrt(R(1,1))*randn(m,L);
v3=sqrt(R(1,1))*randn(m,L);
v5=sqrt(R(1,1))*randn(m,L);
v7=sqrt(R(1,1))*randn(m,L);



 %������������ ������� ������ �����-����� � ��������������� ������� ���������

        xj1=0.25+(width_x-0.25)*rand(1,j1); % samples distribution of the first coordinates 
        xj2=0.25+(width_x-0.25)*rand(1,j1); % samples distribution of the second coordinates
         
%         kk=[atand((width_x-x1)/(length_y-x2)),360-atand(x1/(length_y-x2)),180-atand((width_x-x1)/x2),180+atand(x1/x2)];
        k1=@(xx,yy)atand((width_x-xx)/(length_y-yy));
        k3=@(xx,yy)360-atand(xx/(length_y-yy));
        k5=@(xx,yy)180-atand((width_x-xx)/yy);
        k7=@(xx,yy)180+atand(xx/yy);
        % Range_to_each landmark modeling for every samples
        for j=1:j1

            Sij1(:,j)=getRange_(xj1(1,j),xj2(1,j),k1(xj1(1,j),xj2(1,j)),getWallNum(xj1(1,j),xj2(1,j),k1(xj1(1,j),xj2(1,j)),[width_x;length_y]),[width_x;length_y])+C1d(j);
            
            Sij3(:,j)=getRange_(xj1(1,j),xj2(1,j),k3(xj1(1,j),xj2(1,j)),getWallNum(xj1(1,j),xj2(1,j),k3(xj1(1,j),xj2(1,j)),[width_x;length_y]),[width_x;length_y])+C2d(j);
            
            Sij5(:,j)=getRange_(xj1(1,j),xj2(1,j),k5(xj1(1,j),xj2(1,j)),getWallNum(xj1(1,j),xj2(1,j),k5(xj1(1,j),xj2(1,j)),[width_x;length_y]),[width_x;length_y])+C3d(j);
            
            Sij7(:,j)=getRange_(xj1(1,j),xj2(1,j),k7(xj1(1,j),xj2(1,j)),getWallNum(xj1(1,j),xj2(1,j),k7(xj1(1,j),xj2(1,j)),[width_x;length_y]),[width_x;length_y])+C4d(j);
            
            
%             Sij1(:,j)=sqrt( (xj1(1,j)- (Or1(1)+Er11(j)))^2 + (xj2(1,j)-(Or1(2)+Er12(j)) )^2 + h1^2 ) + C1d(j); 
%             
%             Sij3(:,j)=sqrt( (xj1(1,j)- (Or2(1)+Er21(j)))^2 + (xj2(1,j)-(Or2(2)+Er22(j)) )^2 + h1^2 ) + C2d(j); 
%             
%             Sij5(:,j)=sqrt( (xj1(1,j)- (Or3(1)+Er31(j)))^2 + (xj2(1,j)-(Or3(2)+Er32(j)) )^2 + h1^2 ) + C3d(j); 
%             
%             Sij7(:,j)=sqrt( (xj1(1,j)- (Or4(1)+Er41(j)))^2 + (xj2(1,j)-(Or4(2)+Er42(j)) )^2 + h1^2 ) + C4d(j); 
%         
       
        end 
%
        Sij=[Sij1(1,:);Sij3(1,:);Sij5(1,:);Sij7(1,:)];
        
 %       R1=R1(1:4,1:4);
%         Sij=[Sij1(1,:);Sij2(1,:)];
        xj=[xj1;xj2];
       
posterior=zeros(5,j1);
posterior(1,:)=1/j1;

P=zeros(2,2,4,L);

sigma2=zeros(2,4,L);
figure
for l=1:L %range for every realization


   y1=getRange_(x1(l),x2(l),k1(x1(l),x2(l)),getWallNum(x1(l),x2(l),k1(x1(l),x2(l)),[width_x;length_y]),[width_x;length_y]) + Ci1d(l) + v1(:,l); %������������ ���������
   
   y3=getRange_(x1(l),x2(l),k3(x1(l),x2(l)),getWallNum(x1(l),x2(l),k3(x1(l),x2(l)),[width_x;length_y]),[width_x;length_y]) + Ci2d(l) + v3(:,l); %������������ ���������
   
   y5=getRange_(x1(l),x2(l),k5(x1(l),x2(l)),getWallNum(x1(l),x2(l),k5(x1(l),x2(l)),[width_x;length_y]),[width_x;length_y]) + Ci3d(l) + v5(:,l); %������������ ���������
   
   y7=getRange_(x1(l),x2(l),k7(x1(l),x2(l)),getWallNum(x1(l),x2(l),k7(x1(l),x2(l)),[width_x;length_y]),[width_x;length_y]) + Ci4d(l) + v7(:,l); %������������ ���������
           
       
   
  y=[y1;y3;y5;y7];
 % y=[y3;y4;y5;y6];

 cov__=zeros(2,2,j1);
  for k=1:4 
        
        Exp=0; Exp1=0; Exp2=0; Exp3=0;
        maxLik=0;
        
        %����� �����-�����
       
        for h=1:j1
            
            maxLik(h)=Maximum_Likelihood_calculation_for_LRF(y(k),Sij(k,h),R1(1,1));
            
            cov__(:,:,h)=xj(:,h)*xj(:,h)';
            
            
%              Exp(h)=exp( -0.5*(y-Sij(:,h))'*R1*(y-Sij(:,h)) );
%              
%              Exp21=xj(:,h)*Exp(h);
%                  
%              Exp31=xj(:,h)*xj(:,h)'*Exp(h);
%           
%              Exp1=(1/j1)*Exp(h)+Exp1;
%             
%              Exp2=Exp21+Exp2;
%                        
%              Exp3=Exp31+Exp3;                    
        end
        posterior(k+1,:)=(maxLik.*posterior(k,:))./sum((maxLik.*posterior(k,:)));
        
        X(1,k,l)=sum(xj1.*posterior(k+1,:));
        X(2,k,l)=sum(xj2.*posterior(k+1,:));
        
        for kk=1:j1
                   cov__(:,:,kk)=posterior(k+1,kk).*cov__(:,:,kk);
                                    
        end    
                        
               P(:,:,k,l)=sum(cov__,3)-[X(1,k,l);X(2,k,l)]*[X(1,k,l),X(2,k,l)];
               
               [out]=ResamplingFast(j1,j1,posterior(k+1,:),xj);
               xj(1,:)=out(1,:);
               xj(2,:)=out(2,:);
               xj1=out(1,:);
               xj2=out(2,:);

         
  end 
  plot( xj1,xj2,'.r');
  xlim([0 width_x])
  ylim([0 length_y])
  drawnow();

  sigma2(:,:,l)=(X(:,:,l)-repmat([x(1,l);x(2,l)],1,4));
end

sigma=sqrt(sum(sigma2.^2,3)./(L-1));
P__=sum(P,4)./L;
figure
subplot(211)
plot(linspace(1,m,m),sigma(1,:),'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
subplot(212)
plot(linspace(1,m,m),sigma(2,:),'r',linspace(1,m,m),sqrt(squeeze(P__(2,2,:))),'b');

figure 
hold on
subplot(231)
plot(linspace(1,m,m),abs(sigma(1,:)'-sqrt(squeeze(P__(1,1,:)))),'b');%,'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(234)
plot(linspace(1,m,m),abs(sigma(2,:)'-sqrt(squeeze(P__(2,2,:)))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
