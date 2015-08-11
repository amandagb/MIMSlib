%% Restoration of Metallomic Images 

clear all 
close all 
clc

elem = 'Zn64';
% load('time.csv');
% load([elem '.csv']);

load time
load (elem)
eval(['A = ' elem ';']);


% A = circshift(A, [-133 0]);
mA = mean2(A);
sA = std(A(:));
% A = (A-mA)/sA; %center data

%plot a single line from image
figure
plot(time(9,:),A(9,:))

y = A;


%%

figure
imagesc(A)
title('Raw Image')

%% Set up needed vectors and values

twh = 6;  % t0*(((1/fsh)^(1/3)) - 1); %washout time in no. of pixels
t = (0:1:40)';
fsh = .010; %fraction of signal height defined for washout time
t0 =  twh/(((1/fsh)^(1/3)) - 1);%delay time
s0 =  1;%pi*((.025/2)^2);%1; 
ps = 2;
% tx = 0:1:99;%.1:9.9;
% Ts = tx(11) - tx(1);
s = (s0)*2*(t0^2)./((t0+t).^3);
% [h0,tn] = LA_impulse(twh,fsh,s0,41,1,ps);

h0 = s;
h0 = h0/sum(h0); %normalize 
H = zeros(1,2*length(t)-1);
H(41:81) = h0;
% H = h0';
h0 = H;

M = 1;     % The pre-conditioner "matrix." Set to 1 to avoid figuring that out.

N = length(y);
% C = convmtx(h0,N);    %Convolution matrix - won't be need for 2D.
% C = C(1:N,1:N);
% [Dx Dy] = dgrad([N 1],'same'); % Derivative Op matrix. not needed in 2D.
% d1 = [-1 0 1]; % 1st order approx of 1st derivative
d1 = [-3 4 -1];  % 2nd order approximation of 1st derivative
% d1 = [-1 2 -1];% approx of 2nd derivative
d2 = d1';

alpha = .5;     % initial Regularization parameter
beta = .0001;   % approximation parameter for l_p regularization
gamma = .05;    % exit criteria
p = 1;          % p norm on data fit term if not using Huber
q = 1;          % p norm on Regularization term

x0 =  conv2(y,flipud(fliplr(H)),'same');%ones(size(y));

%CG initial parameters
max_it = 5; 
tol = .7;


%% solve l-p data fit regularization  

err = 1e5; %Set initial error high so we don't exit before starting.
    
[L1,L2] = size(x0);
xhat = x0;

epsF = 1e-5; % Threshold b/w quadratic and linear in Huber norm.
% Compute these once.
fFe = epsF^(p-2);
fFx = ones(L1,L2)*fFe;



tic
    it = 0;
    xhat2 = A;
    it_tot = 0;
    
while err>gamma
    it = it+1 % iteration number
    
% Find weighting matrix Wk at current iteration
    fit = conv2(xhat,H,'same')-y;
    fFx = fit.^2;                   %Huber quadratic section
%     fFx = ones(L1,L2)*fFe;
    locs = abs(fit)>epsF;
    tmp1 = 2*epsF*abs(fit)-epsF^2;  % <-Huber Cost, linear section 
%     tmp1 = abs(fit).^(p-2);       % approx lp;
    fFx (locs) = tmp1(locs);
    tmpF = ((2/p).*fFx);%

% old code to find weighting matrix for Regularization term
%     fRx = zeros(L1,L2);
%     dx = (conv2(xhat,d1,'same')  + conv2(xhat,d2,'same'));
%     locsR = abs(dx)>epsR;
%     tmp2 = abs(dx).^((q-2)/2);
%     fRx (locsR) = tmp2(locsR);
%     tmpR = (2/q).*fRx;%
    
% TV cost
    tmp = (q/2).*( (abs(conv2(xhat,d1,'same')+ conv2(xhat,d2,'same')).^2) + beta ).^(1-q/2);
    

% solve ||y-Cx||_huber + alpha^2 ||Dx||_p with weighted CG
    
    figure
    [xhat2,error, iter, flag] = cgconvtik2(H,tmpF,d1,d2,tmp,xhat,y, alpha,M,max_it, tol);
    
    it_tot = it_tot+iter;
    err = norm(xhat2-xhat,2)/norm(xhat,2);
    
    imagesc(xhat2)  %final result of this iteration
    title(['error = ' num2str(err)])

    disp(['Percent change since last iteration: ',num2str(err*100)])

    % Update estimate
    xhat = xhat2;
    alpha = alpha/2;    % Increase or decrease regularization param
    tol = tol/2;        % decrease CG tolerance
    if(it>5), break;    % Only do 5 iterations
    end
end;
toc

figure
plot(time(9,:),xhat(9,:))


