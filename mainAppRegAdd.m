%% Robust MPC Simple Strategy: Generates the ROA or N-Step Rob Reachable Set Approximate
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
W = Polyhedron('lb',-err_modBound,'ub',err_modBound);                                                         % NET tilde W

%% Form the terminal set and Cost here 
W_Term = Polyhedron('lb',wlb*ones(nx,1),'ub',wub*ones(nx,1));                                          
C = [X.A; zeros(size(U.A,1), nx)]; 
D = [zeros(size(X.A,1), nu); U.A]; 
b = [X.b; U.b]; 
[Xn, Pinf] = term_setRobPar(Anom, Bnom, delAv, delBv, K, X, U, W_Term, Q, R, nx, nu);        

%%% Needed for constraint loop
for i = 1:size(delAv,2)/nx
    setdelA(:,:,i) = delAv(:,(i-1)*nx + 1: i*nx);  
end

for i = 1:size(delBv,2)/nu
    setdelB(:,:,i) = delBv(:,(i-1)*nu + 1: i*nu);  
end

%% Pick this based on what we need 
N_start = 1;            % N or 1. N gives the approx. N-Step robust controllable set. 1 gives the approx ROA.  

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

%% Main Loop Runs start here
%%% vector directions to get inner approximate. Pick anything here. 
dVector{1} = [1;1];
dVector{2} = [0;1];
dVector{3} = [1;0];
dVector{4} = [-1;1];
dVector{5} = [2;-6];
dVector{6} = [2;6];
dVector{7} = [-6;8.2];
dVector{8} = [8.1;-6.2];
dVector{9} = [8.0;-4.439];
%%% check negative directions too
vSign{1} =  1;
vSign{2} = -1;

%% Main Loop 
for Nhor = N_start:N    
    x0feas = []; 
    for i = 1:size(dVector,2)
        for j = 1:2
            [x0feas_out{j}, x0feasNormOut(j)] = FTOCP_add(dVector{i}, vSign{j}, Anom, Bnom, Nhor, X, U, Xn, setdelA, setdelB, W, W_Term, nx, nu, ...
                                                                                          dim_t{Nhor}, matF{Nhor}, matG{Nhor}, matH{Nhor}, mat_c{Nhor}); 
            if x0feasNormOut(j) ~= -inf
                x0feas = [x0feas, x0feas_out{j}];               % add only if feasible 
            end          
        end
    end
    
    x0feasM{Nhor} = x0feas;                                     % storing all feasible points horizon-wise
end

%% Plot the Approx ROA
% It's the union for all Nhors; not the CVX hull of the union! (Although the result in the paper won't vary)
figure; 
for i = 1: length(x0feasM)
    NRC = Polyhedron(x0feasM{i}');                              % corresponding approx. N-step robust controllable set
    plot(NRC, 'color', 'b'); hold on; 
end
