%% Defining System and Constraint Matrices with this function 
%  Monimoy Bujarbaruah

function [Anom,Bnom, epsA, epsB, delAv, delBv, K, A, B, X, U, Xub, Uub, nx, nu, wub,wlb, x_0, Q, R, simsteps, N] = sys_loadNew()

    %%%% Considering two states and one scalar input 
    Anom = [1, 0.15; 0.1, 1];
    Bnom = [0.1; 1.1];                                                          
    nx = size(Anom,2); nu = size(Bnom,2); 

    %%%%%%%%% Get the error vertices matrices
    % these are the \Delta^{(i)}_{A,B} matrices 
    delAv = [ [0, 0.1; 0.1, 0], [0, 0.1; -0.1, 0],...
               [0, -0.1; 0.1, 0], [0, -0.1; -0.1, 0]];                                                                          

    delBv = [[0; -0.1], [0; 0.1], [0.1; 0], [-0.1; 0] ];

    %%%% Fix the error infinity norm limits on the matrices 
    epsA = 0.1; 
    epsB = 0.1; 
    
    %%%% Set the true A and B matrices (satisfy the above bounds)
    A = [1, 0.1; 
           0, 1.0 ]; 
    B = [0; 
         1.0];                    

    %%%%% Weights %%%%%%%%%%%%%%%%%%%%%%%
    Q =  10*eye(nx);
    R =   2*eye(nu);

    %%%%%%%%% choose the feedback gain K here%%%
    %%% NOTICE THAT WE DO (A-BK) LATER
    K = place(Anom, Bnom, [0.72; 0.75]);                      

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simsteps = 100;           % total task length 
    N= 5;                     % try horizon from N to 1

    %%%% Considering constraints of the form -a<=x(i)<=a and -ulb<=u<=uub
    Xlb = -[8; 8];
    Xub = -Xlb; 
    Ulb =  -4; 
    Uub = -Ulb; 

    X = Polyhedron('lb',Xlb,'ub',Xub);
    U = Polyhedron('lb',Ulb,'ub',Uub);            % State and input constraints  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Defining Noise  Bounds 
    wub = 0.1;                                    % Upper bound of additive noise value
    wlb = -0.1;                                   % Lower bound of additive noise value

    %% Starting Condition
    x_0 = [2; -2]; 

end
