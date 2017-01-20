
% 5 plot 200 samples from y

% To plot 200 samples from y, first we should get two handred samples from
% x1 and x2. To do that we should get samples from the marignal
% distribution that we have got in problem 1.
% To sample from a custom distribution the following steps were taken
% 1- get a probability of a group od number from -1 to 1 (interval of X)
% 2- get cdf for this group of numbers
% 3- using interp1 matlab function to get 200 samples from X according to
% the otained cdf
% 4- using y equation in 3 to get 200 samples from y

mypdf=@(x)2/pi*sqrt(1-x.^2); % pdf for x1 and x2 

x=linspace(-1,1,20000); % get group of number in the interval of -1 and 1 
x=mypdf(x); % getting the probability of this numbers
x_norm=x./sum(x); % normalization

figure(1)
plot(linspace(-1,1,20000),x);

cdf_=cumsum(x_norm); % get cdf
figure
plot(linspace(-1,1,20000),cdf_); % check that it starts from zero to 1

X=rand(2,200);

[cdf__, mask] = unique(cdf_);
xq=linspace(-1,1,20000);
xq = xq(mask);

projection = interp1(cdf__, xq, X); % samples according to cdf

figure(1)
hold on
histogram(projection(1,:),'normalization','pdf');


%%from problem 3
% y=[1  -1  X + [2
%    0  2]       3]

y=zeros(2,200);
for i=1:200
y(:,i)=[1 -1;0 2]*projection(:,i)+[2;3];

end

% we can check if the mean and the variance of y is equivalent to our
% analytical solution in problem 3 as follows

mean_y=[mean(y(1,:));mean(y(2,:))];
cov_y=cov(y(1,:),y(2,:));

figure
plot(y(1,:),y(2,:),'.r');
axis equal