function [Xj2,Xj3,B]=ResamplingFast(j1, j2, q, Xj1)

%  clear all;close all; clc;
%  
%  randn('state',0);
% % 
%  j1=3; j2=10000;
% % 
%  Xj1 =[0.0 0.5   1];
%  Exp5=[0.3 0.2 0.5];
%  
%  q=Exp5;

% X=6;b0=1;r=0.1;
% 
% Xj1=randn(1,j1)*b0+X;
% 
% y=5.5+randn*r;
% 
% for i=1:j1 q1(i)=exp( -(0.5/r^2)*(y-Xj1(i))'*(y-Xj1(i)) ); end
% 
% q=q1/sum(q1); 

c1=zeros(j1,1); u = zeros(j2,1);

for i=2:(j1+1);
    c1(i)=c1(i-1)+q(i-1);
end

c=c1(2:end);

Size=size(Xj1);

Xj2=zeros(Size(1),j2); J1=1/j2;

%sum(Exp5)

U=0 + (J1-0).*rand;

k=1;

        for i=1:j2
            
            u(i) = U + J1*(i-1);
            
            while u (i) > c(k)
                k=k+1;
            end
            
                  
            Xj2(:,i)=Xj1(:,k);
            
        end
        
                   
%             for j=0:(j2-1)
%                 if j==0 
%                     Q1=0; Q2=q(1);
%                 else
%                     Q1=Q1+q(j); Q2=Q2+q(j+1);
%                 end
% 
%                 if U(i)>Q1 && U(i)<=Q2
%                    Xj2(i)=Xj1(j+1);
%                    break;
%                    %Exp6(i)=Q2-Q1;
%                 end
% 
%             end
%         end

% x=(-3*b0+X):0.01:(3*b0+X);
% 
% [l,m]=size(x);
% 
% for i=1:m AP1(i)=exp( -(0.5/r^2)*(y-x(i))'*(y-x(i)) - 0.5*( (x(i)-X)^2/(b0^2) )); end
% 
% AP=AP1/sum(AP1);
      
        
% [Xj3,a]=sort([ Xj2 ]);
% 
% j=1;f=1;B=1;
% 
% while f<j1
%     
%     if Xj3(j)==Xj3(j+1)
%         
%         Xj3 (j+1)=[];
%       %  Exp7(j+1)=[];
%         B(j)=B(j)+1;
%         
%     else
%         j=j+1;
%         B=[B 1];
%     end
%        
%     f=f+1;
%     
% end
% 
% Xj3;
% figure(1);hist(Xj1,50);grid;xlim([-3*b0+X 3*b0+X])
% figure(2);hist(Xj2,50);grid;xlim([-3*b0+X 3*b0+X])
% figure(3);stem(Xj3,B);grid;axis([-3*b0+X 3*b0+X 0 max(B)])%'LineWidth',5
% figure(4);plot(x,AP);grid
% figure(5);stem(Xj1,q);grid;xlim([-3*b0+X 3*b0+X])
