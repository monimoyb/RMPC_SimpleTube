%% Finite time optimal control problem: This one used in the time computation part
% Monimoy 


function [feas_flag, cost_flag, v_hor, sol_time] = FTOCP_addTime(x_0, Q, R, Pinf, Anom, Bnom, N, X, U, Xn, setdelA, setdelB, W, W_Term, nx, nu, dim_t, matF, matG, matH, mat_c)

  options = sdpsettings('solver','gurobi','verbose',0);
 
  %% solve the opt problem  
  constraints = []; 
  
   %% Creating Open Loop Optimization Variables for MPC 
    
    M = sdpvar(nu*N,nx*N);        
        for j=1:nu*N
            for k = 2*j-1:nx*N
                    M(j,k) =0;
            end
        end        
     v = sdpvar(nu*N,1);    
     x0feas = x_0;
    
     %% Build the cost 
    x_pred = sdpvar(nx,N+1);
    x_pred(:,1) = x0feas; 
    
    cost_state = (x_pred(:,1))'*Q*(x_pred(:,1)); 
    
    for k=1:N
        x_pred(:,k+1)  = Anom*x_pred(:,k) + Bnom*v(1+(k-1)*nu:k*nu,1);
        if k~=N
        cost_state = cost_state + (x_pred(:,k))'*Q*(x_pred(:,k));
        end
    end
    cost_state = cost_state + (x_pred(:,N+1))'*Pinf*(x_pred(:,N+1));
    
    obj_ol = v'*kron(R,N)*v + cost_state; 
     
   %% Solving the Optimization Problem 
   %%% Cleanly separate the cases here. %%%%%%
   
   if N ==1 
        boldAbar = kron(eye(N),Anom);
        boldBbar = kron(eye(N),Bnom); 
        Fx = blkdiag(kron(eye(N-1), X.A), Xn.A); 
        fx = [kron(ones(N-1,1),X.b); Xn.b]; 
        boldHw = kron(eye(N), W_Term.A);
        boldhw = kron(ones(N,1),W_Term.b); 
        boldHu = kron(eye(N), U.A);
        boldhu = kron(ones(N,1), U.b); 
        %%%variables 
        Lambda = sdpvar(size(Fx,1), size(boldHw,1), 'full'); 
        constraints = [constraints; Lambda>=0];
        gamma = sdpvar(size(boldHw,1), size(boldHu,1), 'full'); 
        constraints = [constraints; gamma >=0];
        constraints = [constraints; gamma'*boldhw <= boldhu - boldHu*v];      
        constraints = [constraints; (boldHu*M)' == boldHw'*gamma];                  %   input con
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% enumerate the set of vertices here 
        for ii = 1:size(setdelA,3)
            bolddelA = kron(eye(N),setdelA(:,:,ii));
            for jj = 1:size(setdelB,3) 
                bolddelB = kron(eye(N),setdelB(:,:,jj));
                constraints = [constraints; Fx*((boldAbar+bolddelA)*x0feas + (boldBbar+bolddelB)*v) + Lambda*boldhw...
                                                                                      <= fx];
            end
        end
   
        constraints = [constraints; Lambda*boldHw == Fx];                               % state con
   
   else      
        %%  Open Loop W Variables  
        Hs=[]; hs =[];   
            for k = 1:N
                polS = W; 
                Hs_ol = polS.A;  hs_ol = polS.b; 
                Hs = blkdiag(Hs, Hs_ol); hs = [hs; hs_ol];
            end
        dim_a = size(Hs,1);                                                  
      
        Z = sdpvar(dim_a,dim_t);  
        constraints = [constraints; matF*v + Z'*hs <= mat_c + matH*x0feas];
        constraints = [constraints; Z>=0];
        constraints = [constraints, matF*M + matG == Z'*Hs];                           % joint con
   end
   
   %% Optimize 
    diagn=solvesdp(constraints, obj_ol, options);
    
    sol_time = diagn.solvertime; 
    feas_flag = diagn.problem; 
   
    
    if feas_flag ~=0
        cost_flag = inf;                             % store high cost if any issue
        v_hor = zeros(nu*N,1);                 % just store some dummy 0's if infeasible anyway
    else
        cost_flag = double(obj_ol);            % store right cost if feasible  
        v_hor = double(v); 
    end
      
   
end