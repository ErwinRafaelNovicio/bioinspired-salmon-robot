function d = solveFEM(K, F)
    n = size(K,1);
    K = sparse(K);

    % quick checks
    if ~issymmetric(K)
        K = (K + K')/2;
    end

    % If RHS is sparse single DOF, convert to sparse
    F = sparse(F);

    % try direct solve (fast & accurate if BC ok)
    try
        d1 = K \ F;
        d=full(d1);
        if norm(K*d1 - F)/max(1,norm(F)) < 1e-8
           
            return
        end
    catch
        % fall through to iterative
    end

    % Setup incomplete cholesky preconditioner
    try
        opts.type = 'ict';
        opts.droptol = 1e-3;
        L = ichol(K, opts);
       
    catch
        L = [];
    end

    % try pcg (K should be SPD-ish)
    tol = 1e-10;
    maxit = 2000;
    try
        [d2, flag, relres] = pcg(K, F, tol, maxit, L, L);
        d=full(d2);
        if flag == 0
            return
           
        end
    catch
        % pcg failed, fallback
    end

    % fallback: small Tikhonov regularization
    alpha = 1e-9 * normest(K);
    d = (K + alpha*speye(n)) \ F;
    d=full(d);

end
