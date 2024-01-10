function [A,n,breaks,dim,jj,C] = splinefit_prep(varargin)
%A butchered version of the splinefit functions:

%   Author: Jonas Lundgren <splinefit@gmail.com> 2010

%basically cut it in half and put everything possible in splinefit_prep
%while everything that depended on data went in splinefit_thin. This was in
%an effort to make things super fast for optimisation code.


% Check number of arguments
error(nargchk(3,7,nargin));

% Check arguments
[x,y,dim,breaks,n,periodic,beta,constr] = arguments(varargin{:});

% Evaluate B-splines
base = splinebase(breaks,n);
pieces = base.pieces;
A = ppval(base,x);

% Bin data
[junk,ibin] = histc(x,[-inf,breaks(2:end-1),inf]); %#ok

% Sparse system matrix
mx = numel(x);
ii = [ibin; ones(n-1,mx)];
ii = cumsum(ii,1);
jj = repmat(1:mx,n,1);
if periodic
    ii = mod(ii-1,pieces) + 1;
    A = sparse(ii,jj,A,pieces,mx);
else
    A = sparse(ii,jj,A,pieces+n-1,mx);
end

% Don't use the sparse solver for small problems
%problems are always small for me!!!
    A = full(A);



    jj = mod(0:pieces+n-2,pieces) + 1;
    
% Compute polynomial coefficients
ii = [repmat(1:pieces,1,n); ones(n-1,n*pieces)];
ii = cumsum(ii,1);
jj2 = repmat(1:n*pieces,n,1);
C = sparse(ii,jj2,base.coefs,pieces+n-1,n*pieces);



%--------------------------------------------------------------------------
function [x,y,dim,breaks,n,periodic,beta,constr] = arguments(varargin)
%ARGUMENTS Lengthy input checking
%   x           Noisy data x-locations (1 x mx)
%   y           Noisy data y-values (prod(dim) x mx)
%   dim         Leading dimensions of y
%   breaks      Breaks (1 x (pieces+1))
%   n           Spline order
%   periodic    True if periodic boundary conditions
%   beta        Robust fitting parameter, no robust fitting if beta = 0
%   constr      Constraint structure
%   constr.xc   x-locations (1 x nx)
%   constr.yc   y-values (prod(dim) x nx)
%   constr.cc   Coefficients (?? x nx)

% Reshape x-data
x = varargin{1};
mx = numel(x);
x = reshape(x,1,mx);

% Remove trailing singleton dimensions from y
y = varargin{2};
dim = size(y);
while numel(dim) > 1 && dim(end) == 1
    dim(end) = [];
end
my = dim(end);

% Leading dimensions of y
if numel(dim) > 1
    dim(end) = [];
else
    dim = 1;
end

% Reshape y-data
pdim = prod(dim);
y = reshape(y,pdim,my);

% Check data size
if mx ~= my
    mess = 'Last dimension of array y must equal length of vector x.';
    error('arguments:datasize',mess)
end

% Treat NaNs in x-data
inan = find(isnan(x));
if ~isempty(inan)
    x(inan) = [];
    y(:,inan) = [];
    mess = 'All data points with NaN as x-location will be ignored.';
    warning('arguments:nanx',mess)
end

% Treat NaNs in y-data
inan = find(any(isnan(y),1));
if ~isempty(inan)
    x(inan) = [];
    y(:,inan) = [];
    mess = 'All data points with NaN in their y-value will be ignored.';
    warning('arguments:nany',mess)
end

% Check number of data points
mx = numel(x);
if mx == 0
    error('arguments:nodata','There must be at least one data point.')
end

% Sort data
if any(diff(x) < 0)
    [x,isort] = sort(x);
    y = y(:,isort);
end

% Breaks
if isscalar(varargin{3})
    % Number of pieces
    p = varargin{3};
    if ~isreal(p) || ~isfinite(p) || p < 1 || fix(p) < p
        mess = 'Argument #3 must be a vector or a positive integer.';
        error('arguments:breaks1',mess)
    end
    if x(1) < x(end)
        % Interpolate breaks linearly from x-data
        dx = diff(x);
        ibreaks = linspace(1,mx,p+1);
        [junk,ibin] = histc(ibreaks,[0,2:mx-1,mx+1]); %#ok
        breaks = x(ibin) + dx(ibin).*(ibreaks-ibin);
    else
        breaks = x(1) + linspace(0,1,p+1);
    end
else
    % Vector of breaks
    breaks = reshape(varargin{3},1,[]);
    if isempty(breaks) || min(breaks) == max(breaks)
        mess = 'At least two unique breaks are required.';
        error('arguments:breaks2',mess);
    end
end

% Unique breaks
if any(diff(breaks) <= 0)
    breaks = unique(breaks);
end

% Optional input defaults
n = 4;                      % Cubic splines
periodic = false;           % No periodic boundaries
robust = false;             % No robust fitting scheme
beta = 0.5;                 % Robust fitting parameter
constr = [];                % No constraints

% Loop over optional arguments
for k = 4:nargin
    a = varargin{k};
    if ischar(a) && isscalar(a) && lower(a) == 'p'
        % Periodic conditions
        periodic = true;
    elseif ischar(a) && isscalar(a) && lower(a) == 'r'
        % Robust fitting scheme
        robust = true;
    elseif isreal(a) && isscalar(a) && isfinite(a) && a > 0 && a < 1
        % Robust fitting parameter
        beta = a;
        robust = true;
    elseif isreal(a) && isscalar(a) && isfinite(a) && a > 0 && fix(a) == a
        % Spline order
        n = a;
    elseif isstruct(a) && isscalar(a)
        % Constraint structure
        constr = a;
    else
        error('arguments:nonsense','Failed to interpret argument #%d.',k)
    end
end

% No robust fitting
if ~robust
    beta = 0;
end

% Check exterior data
h = diff(breaks);
xlim1 = breaks(1) - 0.01*h(1);
xlim2 = breaks(end) + 0.01*h(end);
if x(1) < xlim1 || x(end) > xlim2
    if periodic
        % Move data inside domain
        P = breaks(end) - breaks(1);
        x = mod(x-breaks(1),P) + breaks(1);
        % Sort
        [x,isort] = sort(x);
        y = y(:,isort);
    else
        mess = 'Some data points are outside the spline domain.';
        warning('arguments:exteriordata',mess)
    end
end

% Return
if isempty(constr)
    return
end

% Unpack constraints
xc = [];
yc = [];
cc = [];
names = fieldnames(constr);
for k = 1:numel(names)
    switch names{k}
        case {'xc'}
            xc = constr.xc;
        case {'yc'}
            yc = constr.yc;
        case {'cc'}
            cc = constr.cc;
        otherwise
            mess = 'Unknown field ''%s'' in constraint structure.';
            warning('arguments:unknownfield',mess,names{k})
    end
end

% Check xc
if isempty(xc)
    mess = 'Constraints contains no x-locations.';
    error('arguments:emptyxc',mess)
else
    nx = numel(xc);
    xc = reshape(xc,1,nx);
end

% Check yc
if isempty(yc)
    % Zero array
    yc = zeros(pdim,nx);
elseif numel(yc) == 1
    % Constant array
    yc = zeros(pdim,nx) + yc;
elseif numel(yc) ~= pdim*nx
    % Malformed array
    error('arguments:ycsize','Cannot reshape yc to size %dx%d.',pdim,nx)
else
    % Reshape array
    yc = reshape(yc,pdim,nx);
end

% Check cc
if isempty(cc)
    cc = ones(size(xc));
elseif numel(size(cc)) ~= 2
    error('arguments:ccsize1','Constraint coefficients cc must be 2D.')
elseif size(cc,2) ~= nx
    mess = 'Last dimension of cc must equal length of xc.';
    error('arguments:ccsize2',mess)
end

% Check high order derivatives
if size(cc,1) >= n
    if any(any(cc(n:end,:)))
        mess = 'Constraints involve derivatives of order %d or larger.';
        error('arguments:difforder',mess,n-1)
    end
    cc = cc(1:n-1,:);
end

% Check exterior constraints
if min(xc) < xlim1 || max(xc) > xlim2
    if periodic
        % Move constraints inside domain
        P = breaks(end) - breaks(1);
        xc = mod(xc-breaks(1),P) + breaks(1);
    else
        mess = 'Some constraints are outside the spline domain.';
        warning('arguments:exteriorconstr',mess)
    end
end

% Pack constraints
constr = struct('xc',xc,'yc',yc,'cc',cc);


%--------------------------------------------------------------------------
function pp = splinebase(breaks,n)
%SPLINEBASE Generate B-spline base PP of order N for breaks BREAKS

breaks = breaks(:);     % Breaks
breaks0 = breaks';      % Initial breaks
h = diff(breaks);       % Spacing
pieces = numel(h);      % Number of pieces
deg = n - 1;            % Polynomial degree

% Extend breaks periodically
if deg > 0
    if deg <= pieces
        hcopy = h;
    else
        hcopy = repmat(h,ceil(deg/pieces),1);
    end
    % to the left
    hl = hcopy(end:-1:end-deg+1);
    bl = breaks(1) - cumsum(hl);
    % and to the right
    hr = hcopy(1:deg);
    br = breaks(end) + cumsum(hr);
    % Add breaks
    breaks = [bl(deg:-1:1); breaks; br];
    h = diff(breaks);
    pieces = numel(h);
end

% Initiate polynomial coefficients
coefs = zeros(n*pieces,n);
coefs(1:n:end,1) = 1;

% Expand h
ii = [1:pieces; ones(deg,pieces)];
ii = cumsum(ii,1);
ii = min(ii,pieces);
H = h(ii(:));

% Recursive generation of B-splines
for k = 2:n
    % Antiderivatives of splines
    for j = 1:k-1
        coefs(:,j) = coefs(:,j).*H/(k-j);
    end
    Q = sum(coefs,2);
    Q = reshape(Q,n,pieces);
    Q = cumsum(Q,1);
    c0 = [zeros(1,pieces); Q(1:deg,:)];
    coefs(:,k) = c0(:);
    % Normalize antiderivatives by max value
    fmax = repmat(Q(n,:),n,1);
    fmax = fmax(:);
    for j = 1:k
        coefs(:,j) = coefs(:,j)./fmax;
    end
    % Diff of adjacent antiderivatives
    coefs(1:end-deg,1:k) = coefs(1:end-deg,1:k) - coefs(n:end,1:k);
    coefs(1:n:end,k) = 0;
end

% Scale coefficients
scale = ones(size(H));
for k = 1:n-1
    scale = scale./H;
    coefs(:,n-k) = scale.*coefs(:,n-k);
end

% Reduce number of pieces
pieces = pieces - 2*deg;

% Sort coefficients by interval number
ii = [n*(1:pieces); deg*ones(deg,pieces)];
ii = cumsum(ii,1);
coefs = coefs(ii(:),:);

% Make piecewise polynomial
pp = mkpp(breaks0,coefs,n);


%--------------------------------------------------------------------------
function B = evalcon(base,constr,periodic)
%EVALCON Evaluate linear constraints

% Unpack structures
breaks = base.breaks;
pieces = base.pieces;
n = base.order;
xc = constr.xc;
cc = constr.cc;

% Bin data
[junk,ibin] = histc(xc,[-inf,breaks(2:end-1),inf]); %#ok

% Evaluate constraints
nx = numel(xc);
B0 = zeros(n,nx);
for k = 1:size(cc,1)
    if any(cc(k,:))
        B0 = B0 + repmat(cc(k,:),n,1).*ppval(base,xc);
    end
    % Differentiate base
    coefs = base.coefs(:,1:n-k);
    for j = 1:n-k-1
        coefs(:,j) = (n-k-j+1)*coefs(:,j);
    end
    base.coefs = coefs;
    base.order = n-k;
end

% Sparse output
ii = [ibin; ones(n-1,nx)];
ii = cumsum(ii,1);
jj = repmat(1:nx,n,1);
if periodic
    ii = mod(ii-1,pieces) + 1;
    B = sparse(ii,jj,B0,pieces,nx);
else
    B = sparse(ii,jj,B0,pieces+n-1,nx);
end


%--------------------------------------------------------------------------
function [Z,u0] = solvecon(B,constr)
%SOLVECON Find a particular solution u0 and null space Z (Z*B = 0)
%         for constraint equation u*B = yc.

yc = constr.yc;
tol = 1000*eps;

% Remove blank rows
ii = any(B,2);
B2 = full(B(ii,:));

% Null space of B2
if isempty(B2)
    Z2 = [];
else
    % QR decomposition with column permutation
    [Q,R,dummy] = qr(B2); %#ok
    R = abs(R);
    jj = all(R < R(1)*tol, 2);
    Z2 = Q(:,jj)';
end

% Sizes
[m,ncon] = size(B);
m2 = size(B2,1);
nz = size(Z2,1);

% Sparse null space of B
Z = sparse(nz+1:nz+m-m2,find(~ii),1,nz+m-m2,m);
Z(1:nz,ii) = Z2;

% Warning rank deficient
if nz + ncon > m2
	mess = 'Rank deficient constraints, rank = %d.';
	warning('solvecon:deficient',mess,m2-nz);
end

% Particular solution
u0 = zeros(size(yc,1),m);
if any(yc(:))
    % Non-homogeneous case
	u0(:,ii) = yc/B2;
    % Check solution
	if norm(u0*B - yc,'fro') > norm(yc,'fro')*tol
        mess = 'Inconsistent constraints. No solution within tolerance.';
        error('solvecon:inconsistent',mess)
	end
end


%--------------------------------------------------------------------------
function u = lsqsolve(A,y,beta)
%LSQSOLVE Solve Min norm(u*A-y)

% Avoid sparse-complex limitations
if issparse(A) && ~isreal(y)
    A = full(A);
end

% Solution
u = y/A;

% Robust fitting
if beta > 0
    [m,n] = size(y);
    alpha = 0.5*beta/(1-beta)/m;
    for k = 1:3
        % Residual
        r = u*A - y;
        rr = r.*conj(r);
        rrmean = sum(rr,2)/n;
        rrmean(~rrmean) = 1;
        rrhat = (alpha./rrmean)'*rr;
        % Weights
        w = exp(-rrhat);
        spw = spdiags(w',0,n,n);
        % Solve weighted problem
        u = (y*spw)/(A*spw);
    end
end

