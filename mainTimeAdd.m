%% Robust MPC Simple Strategy: Generates the Times 
% This method lumps up model mismatch component as additive uncertainty for N > 1
% BUT, using exact uncertainty for N=1 and the terminal set 
% Author: Monimoy Bujarbaruah

%%
clear all
close all
clc
yalmip 'clear'

%% MPC Controller Parameters
[Anom,Bnom, epsA, epsB, delAv, delBv, K, A, B, X, U, Xub, Uub, nx, nu, wub,wlb, x_0, Q, R, simsteps, N] = sys_loadNew();

%% Form the net additive error bound here
err_modBound = epsA*Xub + epsB*Uub + wub;   
W = Polyhedron('lb',-err_modBound,'ub',err_modBound);                                                     % NET tilde W

%% Form the terminal set and Cost here 
W_Term = Polyhedron('lb',wlb*ones(nx,1),'ub',wub*ones(nx,1));                                          
C = [X.A; zeros(size(U.A,1), nx)]; 
D = [zeros(size(X.A,1), nu); U.A]; 
b = [X.b; U.b]; 
[Xn, Pinf] = term_setRobPar(Anom, Bnom, delAv, delBv, K, X, U, W_Term, Q, R, nx, nu);        

%% For closed loop trajectory recording  
x_cl = zeros(nx,simsteps+1);
u_cl = zeros(nu,simsteps); 

%% Needed for constraint loop
for i = 1:size(delAv,2)/nx
    setdelA(:,:,i) = delAv(:,(i-1)*nx + 1: i*nx);  
end

for i = 1:size(delBv,2)/nu
    setdelB(:,:,i) = delBv(:,(i-1)*nu + 1: i*nu);  
end

%% N problems to be solved here at any step
N_start = 1;             

%% Form all the required stuff for all horizon options
for Nhor = N_start:N
    dim_ttmp = size(C,1)*Nhor + size(Xn.A,1);  
    [capA, capE, capB, capC, capD, Aw_batch, Bu_batch, A_batch] = obtain_matR(Anom, Bnom, C, D, Xn, nx, nu, Nhor, dim_ttmp);
    matF{Nhor} = capC*capB + capD; 
    matG{Nhor} = capC*capE; 
    matH{Nhor} = -capC*capA; 
    mat_c{Nhor} = [kron(ones(Nhor,1),b); Xn.b];   
    dim_t{Nhor} = dim_ttmp;
end

%% Getting Started to Solve the MPC Problems and Record Time
x_cl(:,1) = x_0;                                                  
N_start = 1;

for i = 1:simsteps                                                                         % will record these times in closed-loop
     x_init = x_cl(:,i);    
     
     for Nhor = N_start:N
            [feas_flag(Nhor), cost_flag(Nhor), v_horN{Nhor}, sol_time(Nhor, i)] = FTOCP_addTime(x_init, Q, R, Pinf, Anom, Bnom, Nhor, X, U, Xn, setdelA, setdelB, W, W_Term, nx, nu, ...
                                                                                             dim_t{Nhor}, matF{Nhor}, matG{Nhor}, matH{Nhor}, mat_c{Nhor}); 
     end
        
    % pick best cost 
    [value_mincost, ind_mincost] = min(cost_flag); 
    
    if value_mincost == inf
        disp('*** infeasible problem***')                                                  % should not happen if no bug
        break
    else
         v_hor = v_horN{ind_mincost};
    end
     
  %% Obtaining Closed Loop Parameters 
    u_cl(:,i) = v_hor(1:nu,:);                                                             % Closed loop control 
  %% Actual System Simulation in Closed Loop 
    w = wlb + (wub-wlb)*rand(nx,1);                                               
    x_cl(:,i+1) = A*x_cl(:,i) + B*u_cl(:,i) +  w;      
    yalmip 'clear'   
    i
end


%% Get the mean solver time 
online_time = mean(sol_time,2);                                                            % For each horizon
