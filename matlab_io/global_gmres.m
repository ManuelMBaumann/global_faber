function [X, error, iter, resvec,flag] = global_gmres( A, X, C, M, restrt, max_it, tol )

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag] = gmres( A, x, b, M, restrt, max_it, tol )
%
% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         restrt   INTEGER number of iterations between restarts
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

    iter = 0;                                         % initialization
    flag = 0;
    b = 1;
    [n,b] = size(C);
    cnrm2 = norm( C,'fro' );

    if  ( cnrm2 == 0.0 ), cnrm2 = 1.0; end

    r = M \ ( C-(A*X+X*A') );
    error = norm( r ,'fro') / cnrm2;
    resvec = error;
    if ( error < tol ) return, end

                              % initialize workspace
    m = restrt;
    V = zeros(n,(m+1)*b);
    H = zeros(m+1,m);
    cs = zeros(m+1,1);
    sn = zeros(m+1,1);
    e1    = zeros(m+1,1);
    e1(1) = 1.0;
    
    for iter = 1:max_it,                              % begin iteration
    
        V(:,1:b) = r / norm( r );
        s = norm( r )*e1;
        for i = 1:m,                                   % construct orthonormal
            w = M \ (A*V(:,(i-1)*b+1:i*b) + V(:,(i-1)*b+1:i*b)*A');                         % basis using Gram-Schmidt
            for k = 1:i,
                H(k,i)= fro_dot(w,V(:,(k-1)*b+1:k*b));  %w'*V(:,k);
                w = w - H(k,i)*V(:,(k-1)*b+1:k*b);
            end
            H(i+1,i) = norm( w, 'fro');
            V(:,i*b+1:(i+1)*b) = w / H(i+1,i);
            for k = 1:i-1,                              % apply Givens rotation
                temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
                H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
                H(k,i)   = temp;
            end
            [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
            temp   = cs(i)*s(i);                        % approximate residual norm
            s(i+1) = -sn(i)*s(i);
            s(i)   = temp;
            H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
            H(i+1,i) = 0.0;
            error  = abs(s(i+1)) / cnrm2
            resvec = [resvec;error];
            if ( error <= tol ),                        % update approximation
                y = H(1:i,1:i) \ s(1:i);                 % and exit
                    %X = X + V(:,1:i)*y;
                X = X+start_product(V,y,b,i);
                break;
            end
        end

    if ( error <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    %X = C + V(:,1:m)*y;                            % update approximation
    X = X+start_product(V,y,b,m);
    r = M \ ( C-(A*X+X*A') );                              % compute residual
    s(i+1) = norm(r);
    error = s(i+1) / cnrm2;                        % check convergence
    resvec = [resvec;error];
    if ( error <= tol ), break, end;
    end

    if ( error > tol ) flag = 1; end;                 % converged

    % END of gmres.m
