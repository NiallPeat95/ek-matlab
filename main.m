%%
clc
diary('ek_base_output.out');
%% ECON 597: Empirical Methods
% 
% Empirical Methods Project - Spring 2011
% Revisiting computational issues of Eaton and Kortum (2002)
% 
% By Bernardo Diaz de Astarloa

% Main file:
% - Run EK iterative procedure with Newton-Raphson step
% - Run alternative algorithms in the Newton-Raphson step
% - KNITRO implementation

% - This file provides baseline results for counterfactuals

tic

%% Arrange data

clear all
run getdata

%% Parameters 

clear all
clc

% Global variables are only parameters which do not change
global a b theta t d n L Y g 

load DATA

n = 19;                              % Number of countries
a = myalpha;                
b  = 0.21221;
theta = 8.28;                        % Shape of Pareto dist
relaw = aw./aw(n,1);                 % Relative wage to US
t = exp( b*(S+theta*(log(relaw))) ); % Technology: see Table VI
d = D;                               % Geographic barriers dni^(-theta)
L = l;                               % Labor
Y = y;                               % Income
elast = 0.1;                         % Elasticity of wage adjustment

% Gamma constant (as in the original code)
g = (b^(-b))*((1-b)^(-(1-b)));                    
g = g^(theta);                                    


%%  Original EK (2002) approach:
%   Solve for the equilibrium using an iterative procedure with a
%   Newton-Raphson intermediate step.

disp('Newton iterative procedure starts here')

w = aw;                     % Set initial wages at data

%%Alternative initial wage
%w = 100*ones(n,1);

mytolw = 1e-9;              
myiter = 20000;
tol = 1;
i = 1;

while (tol>mytolw && i<myiter)
    
    i
    
    % Set bounds (implied by the model) for prices
     pbounds = bounds(w);                   
     p0 = 0.5*(pbounds(:,1)+pbounds(:,2));  % Half way between min and max
  
    %Alternative: fsolve;
    options=optimset('Display','iter','Jacobian','on','TolFun',1e-6);
    [pval,fval,exitflag,jacobian,output] = ...
    fsolve(@(p) prices(p,w),p0,options);
    
    p = pval;                       % Prices
    
    Pi = shares(p,w);               % Get (transpose of) trade shares
    
    yl = wbill(Pi);                 % Get wage bill
    
    exl = (yl./w - l)./l;           % Implied excess labor demand
   
    % Update wages using excess demand for labor
    wnew = w.*(ones(n,1) + exl*elast);      
    w = wnew;
    tol = max(abs(exl)); 
    i = i+1;
    
    [p,w,exl];
    
end

disp('--- Exit flag ---')
if i == myiter;
    disp('Maximum number of iterations reached');
else
    disp('Loop converged');
end

disp('----------------------')
disp('Time consumed: ')
toc
disp('----------------------')


eqw = w;                                % Equilibrium wages
eqrw = eqw./eqw(n,1);                   % Relative to the US

% Measure of "fit" in terms of wages
mse = mean((log(eqw) - log(aw)).^2);
rmse = sqrt(mse);
infnorm = max(abs(eqw-aw));

pbounds = bounds(eqw);
p0 = 0.5*(pbounds(:,1)+pbounds(:,2));
p = newton(@(p) prices(p,eqw),p0);

eqp = p.^(-1/theta);                    % Manufacturing prices
eqrp = eqp./eqp(n,1);                   % Relative to the US

realw = eqw./(eqp.^a);                  % Real wages

basew = [eqw,realw,eqrw,aw];
basep = [eqp,eqrp];

tradeshares = shares(p,eqw)';           % Predicted trade shares

% Total manufacturing demand 
ym = ((1-b)/b)*yl + a*Y;

% Implied exports and imports
trademat = tradeshares.*repmat(ym,1,n);
predxnn = diag(trademat);                 % Predicted domestic spending
predexps = sum(trademat)' - predxnn;      % Predicted exports 
predimps = sum(trademat,2) - predxnn;     % Predicted imports

basetrade = [predexps,predimps];

% Welfare
welfare = Y.*(eqp.^(-a));

disp('----------- Summary of results: Newton -----------')


disp('Model fit (wages)')
disp(['Root mean squared error: ',num2str(rmse)])
disp(['Max of absolute deviations: ',num2str(infnorm)])

disp('Equilibrium abs, real, and relative wages, and data')
disp(basew)

disp('Equilibrium Prices (absolute & relative)') 
disp(basep)

disp('Check trade shares sum to 1')
sum(tradeshares,2)

% disp('Predicted and actual exports')
% disp([exps,predexps])
% 
% disp('Predicted and actual imports')
% disp([imps,predimps])

% Predicted trade shares
figure(1)
mesh(tradeshares); 
figure(gcf)

% Actual trade shares
figure(2)
mesh(Pi_in); 
figure(gcf)

save baseline basew basep basetrade welfare

 toc
diary('off');

