function [x,resvec,restime,spmv_history] = fastgmres(A,b,tol,maxit,M1,M2,x0,opts)
% FASTGMRES   Flexible and sketched GMRES with Arnoldi truncation and SpMV tracking
%
% [x,resvec,restime,spmv_history] = fastgmres(A,b,tol,maxit,M1,M2,x0,opts)

tstart = tic;
restime = []; 
resvec = [];  
spmv_history = []; % Tracks cumulative SpMVs at the end of each outer iteration

if nargin==8 && isfield(opts,'verbose') 
    verbose = opts.verbose;
else
    verbose = 1;
end
if nargin==8 && isfield(opts,'maxtime') 
    maxtime = opts.maxtime;
else
    maxtime = inf;
end
if nargin==8 && isfield(opts,'t') 
    t = opts.t;
else
    t = 0;
end
if nargin==8 && isfield(opts,'Sfun') 
    Sfun = opts.Sfun;                
else
    Sfun = @(N,s) clarkson_woodruff_inner(N,s);
end
if nargin==8 && isfield(opts,'m') 
    m = opts.m;
else
    m = 500;
end
if nargin==8 && isfield(opts,'s') 
    s = opts.s;
else
    s = 2*m;
end
assert(m<s,'embedding dimension s must be > m inner iterations')
if nargin==8 && isfield(opts,'cndtol') 
    cndtol = opts.cndtol;
else
    cndtol = 1e15;
end

% Set up the global SpMV counter
global TOTAL_SPMV_COUNT;
if isempty(TOTAL_SPMV_COUNT)
    TOTAL_SPMV_COUNT = 0;
end

% Wrap operators to automatically increment the global SpMV counter
A_counted = @(x) count_and_eval(A, x);

if isnumeric(A)
    A_counted = @(x) count_and_eval(@(y) A*y, x);
end

N = size(b,1);
if nargin >= 7 && ~isempty(x0)
    r0 = b - A_counted(x0);
else
    x0 = zeros(N,1);
    r0 = b;
end

if nargin >=5 && ~isempty(M1) && ~isempty(M2)
    hA = @(x) M2\(M1\(A_counted(x)));
    r0 = M2\(M1\r0);
    rhs = M2\(M1\b); 
    nrmrhs = norm(rhs); 
else
    hA = A_counted;
    rhs = b;
    nrmrhs = norm(rhs);
end

if nargin < 4
    maxit = 500;
end

W = zeros(N,maxit+1);
Z = zeros(N,maxit);
H = zeros(maxit+1,maxit);
nrmr0 = norm(r0);
W(:,1) = r0/nrmr0;

% Define inner solver using wrapped counted operator
solver = @(rhs) sgmres_inner(A_counted, rhs, m, cndtol, t, s, Sfun, verbose);

for j = 1:maxit
    if verbose >=2, fprintf('  FGMRES iteration = %d\n',j); end

    wvec = W(:,j);
    z = solver(wvec);
    Z(:,j) = z;      
    w = hA(z); % Outer SpMV (+1 SpMV)

    for i = 1:j
        H(i,j) = W(:,i)'*w;
        w = w - W(:,i)*H(i,j);
    end

    H(j+1,j) = norm(w);
    W(:,j+1) = w/H(j+1,j);

    y = H(1:j+1,1:j)\(nrmr0*eye(j+1,1));
    resvec(j) = norm(nrmr0*eye(j+1,1) - H(1:j+1,1:j)*y);
    
    % Record the cumulative SpMV count at the end of this outer iteration
    spmv_history(j) = TOTAL_SPMV_COUNT;

    if verbose >= 2, fprintf('  FGMRES residual  = %5.3e (Cumulative SpMVs: %d)\n\n', resvec(j), spmv_history(j)); end
    restime(j) = toc(tstart);

    if resvec(j)/nrmrhs < 0.99*tol || restime(j) > maxtime
        break
    end
end

x = x0 + Z(:,1:j)*y;
res_true = norm(rhs - hA(x))/nrmrhs;
if res_true > tol
    warning('Target residual NOT reached after %d iterations: %5.3e > %5.3e \n         Consider increasing maxit or restart with x as initial guess.',j,res_true,tol);
end
if verbose && res_true <= tol
    fprintf('fastgmres reached target residual after %d iterations: %5.3e <= %5.3e (Total SpMVs: %d)\n', j, res_true, tol, TOTAL_SPMV_COUNT);
end

end

% Helper function to automatically count matrix evaluations
function y = count_and_eval(A_op, x)
    global TOTAL_SPMV_COUNT;
    TOTAL_SPMV_COUNT = TOTAL_SPMV_COUNT + 1;
    y = A_op(x);
end

function [x,res,sres,cnd,scnd] = sgmres_inner(A, b, m, cndtol, t, s, Sfun, verbose)
N = size(b,1);
hS = Sfun(N,s);

res = norm(b); sres = []; cnd = []; scnd = [];
H = zeros(m+1,m); V = zeros(N,m+1);
v = b/norm(b); V(:,1) = v;
SAV = zeros(s,m+1);
SV = zeros(s,m+1);
SV(:,1) = hS(V(:,1));
Sb = hS(b);
Q = []; R = [];
if verbose >= 2, fprintf('  running sgmres '); end
for j = 1:m
    if verbose >= 2, fprintf('.'); end
    w = A(v); % Inner SpMV (+1 SpMV)
    SAV(:,j) = hS(w);
    [Q,R] = qrupdate_gs_inner(SAV(:,1:j),Q,R);

    if j > 1 && cond(R) > cndtol
        break
    end

    y = R\(Q'*Sb);   
    for i = max(j-t+1,1):j
        H(i,j) = V(:,i)'*w;
        w = w - V(:,i)*H(i,j);
    end
    H(j+1,j) = norm(w);
    v = w/H(j+1,j);
    V(:,j+1) = v;
end
if verbose >= 2, fprintf('\n'); end
x = V(:,1:length(y))*y;
end 

function S = clarkson_woodruff_inner(N, s)
i = randi(s,N,1);
j = 1:N;
el = 2*randi(2,N,1)-3;
S = sparse(i,j,el,s,N,N);
S = @(x) S*x;
end

function [Q1,R1] = qrupdate_gs_inner(A,Q,R)
[m,n] = size(A);
k = size(Q,2);
assert(m>n & k<n, 'A needs to be tall and size(Q,2)<size(A,2)')
Q1 = zeros(m,n);
R1 = zeros(n,n);
Q1(1:m,1:k) = Q;
R1(1:k,1:k) = R;
for j = k+1:n
    w = A(:,j);
    for reo = 0:1
        for i = 1:j-1
            r = Q1(:,i)'*w;
            R1(i,j) = R1(i,j) + r;
            w = w - Q1(:,i)*r;
        end
    end
    R1(j,j) = norm(w);
    Q1(:,j) = w/R1(j,j);
end
end