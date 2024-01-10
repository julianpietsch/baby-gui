function pp = splinefit_thin(A,n,breaks,dim,jj,C,y)
%A butchered version of the splinefit functions:

%   Author: Jonas Lundgren <splinefit@gmail.com> 2010

%basically cut it in half and put everything possible in splinefit_prep
%while everything that depended on data went in splinefit_thin. This was in
%an effort to make things super fast for optimisation code.



% Solve

% Solve Min norm(u*A-y)
u = y/A;


% Periodic expansion of solution
%was an 'if periodic' loop here, but for my purposes the solution is always
%periodic!!!
u = u(:,jj);


% Compute polynomial coefficients

coefs = u*C;
coefs = reshape(coefs,[],n);

% Make piecewise polynomial
pp = mkpp(breaks,coefs,dim);

end

