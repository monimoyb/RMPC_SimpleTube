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
[Anom,Bnom, epsA, epsB, delAv, delBv, K, A, B, X, U, Xlb, Xub, Ulb, Uub, nx, nu, wub,wlb, x_0, Q, R, simsteps, N] = sys_loadNew();
%% Form the net additive error bound here
% ||A||_1 \leq sqrt(m) ||A||_2 for m*n matrix A
err_modBound1 = epsA*Xub + epsB*Uub + wub;  
err_modBound2 = -epsA*Xlb - epsB*Ulb - wlb;                                                                   % can be ASYMMETRIC BOUND
err_modBound = max(err_modBound1, err_modBound2); 
W = Polyhedron('lb',-err_modBound,'ub',err_modBound);                                                     % NET W

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
N_start = 1;            % N or 1. N gives the approx. N-Step robust reachable set. 1 gives the approx ROA.  
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
count_inf = 0;                      % counts the number of infeasible points
xfeas = [];
xinfs = []; 
%%% vector directions to get inner approximate. Pick anything here. 
dVector{1} =  [1;1];
dVector{2} = [ 0;1];
dVector{3} = [ 1;0];
dVector{4} = [-1;1];
dVector{5} = [ 2;-6];
dVector{6} = [2;6];
dVector{7} = [-6;8.2];
dVector{8} = [8.1;-6.2];
dVector{9} = [8.0;-4.439];
vSign{1}    =  1;
vSign{2}    = -1;
x0feas = [];
counter = 1;
idxBestNorm = [];
%%
for i = 1:size(dVector,2)
    for j =1:2
        for Nhor = N_start:N
            [x0feas_out{Nhor}, x0feasNormOut(Nhor)] = FTOCP_add(dVector{i}, vSign{j}, Anom, Bnom, Nhor, X, U, Xn, setdelA, setdelB, W, W_Term, nx, nu, ...
                                                                                              dim_t{Nhor}, matF{Nhor}, matG{Nhor}, matH{Nhor}, mat_c{Nhor}); 
        end
        
        % pick best cost 
        [~, ind_maxNorm] = max(x0feasNormOut); 
        x0feas_normout{i,j} = x0feasNormOut;
        
        idxBestNorm = [idxBestNorm, ind_maxNorm];
    
        x0feas = [x0feas, x0feas_out{ind_maxNorm}];
    
    end
end


%% Plot the ROA (Actually, only the union; not the CVX hull of the union! But plotting the CVX hull here gives the correct set too.)
Xfeas = Polyhedron(x0feas');
figure
plot(Xfeas, 'color', 'b')
