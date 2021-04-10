%% Finite time optimal control problem: This one used in the ROA computation part
% Monimoy Bujarbaruah


function [x0feas_out, x0feasNormOut] = FTOCP_add(dVector, vSign, Anom, Bnom, N, X, U, Xn, setdelA, setdelB, W, W_Term, nx, nu, dim_t, matF, matG, matH, mat_c)

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
     x0feas = sdpvar(nx,1);
     
   %% Solving the Optimization Problem 
   %%% Cleanly separate the cases here. 
   
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
   
        constraints = [constraints; Lambda*boldHw == Fx];                           % state con
   
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
        constraints = [constraints, matF*M + matG == Z'*Hs];                          % joint con (See www.sciencedirect.com/science/article/pii/S0005109806000021)
   end
   
   constraints = [constraints; dVector(2)*x0feas(1) == dVector(1)*x0feas(2)]; 
   constraints = [constraints; X.A*x0feas <= X.b];                                    % common for both cases
   
    cost = vSign*dVector'*x0feas;  
   
    diagn=solvesdp(constraints, cost, options);
   
    feas_flag = diagn.problem; 
   
    
    if feas_flag ~=0
        cost_flag = inf;                    % store high cost if any issue
        x0feas_out = double(x0feas);
        x0feasNormOut = -inf;               % signal of infeasibility!
    else
        cost_flag = double(cost);           % store right cost if feasible  
        v_hor = double(v); 
        x0feas_out = double(x0feas);
        x0feasNormOut = norm(x0feas_out,2);
    end
      
   
end
