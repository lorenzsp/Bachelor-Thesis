% REGRESSIONE_LINEARE returns the fit parameters for a linear fit of the form  y = m*x + b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: REGRESSIONE_LINEARE returns the textbook-formulas based parameters for a linear fit 
% of the form  y = m*x + b.
%
% CALL:
% out = regressione_lineare(x, y, dy, varargin)
% out = regressione_lineare(x, y, [], varargin)
% out = regressione_lineare(x, y)
% out = regressione_lineare(x, y, dy, 'dx', dx)
%
% Where the inputs are:
% x  a vector of data of the independent physical quantity
% y  a vector of data of the dependent physical quantity
% dy a vector of uncertainty of the independent physical quantity
%
% The output is a struct whose fields are
%
% m     the estimated value for the slope of the line
% b     the estimated value for the intercept of the line
% dm    the estimated value for the uncertainty of the slope of the line
% db    the estimated value for the uncertainty of the intercept of the line
% C     the estimated value for the covariance of the line parameters
% chi2  the overall discrepancy between the best straight line and the experimental points
%       NB gives "reduced chi2", normalized per degrees of freedom
% dof   the number of degrees of freedom of the procedure of best fit
%
%
% Notes:
% x and y must have the same number of elements
% dy must have the same number of elements of y, or 1, or none
%
% If errors dy are empty, regressione_lineare gives equal weight to all points and
% estimates the parameters uncertainties based on the assumption that Chi^2 = 1
%
% Optional parameters are:
%
% 'dx'      a vector or single value (assumed equal on all points) of error in x
%           follows. If specified, the fit is first computed without
%           accounting for x errors. Then the slope is used to convert the
%           x error in an y one, dy is updated accordingly and the fit is
%           performed again.
%

function out = regressione_lineare (x, y, dy, varargin)
  
  % Default quantities
  dx = [];
  noerrors    = false;
  
  % Check we have vectors or columns
  if ~isvector(x) || ~isvector(y)
    error('Error: this function only works with vectors or scalars.');
  end
  
  % If the user did not pass the dy quantity, set it to empty
  if ~exist('dy', 'var')
    dy = [];
  end
  
  % Make sure all inputs are columns
  if isrow(x)
    x  =  x';
  end
  if isrow(y)
    y  =  y';
  end
  if isrow(dy)
    dy  =  dy';
  end
  
  if numel(dy) == 1
    dy = dy * ones(size(y));
  end
  
  if isempty(dy)
    noerrors = true;
  end
  
  % Number of points
  N = numel(x);
  
  % Search for optional arguments
  if ~isempty(varargin)
    jj = 1;
    while jj <= length(varargin)
      switch lower(varargin{jj})
        case 'dx'
          dx = cell2mat(varargin(jj+1));
          jj = jj + 1;
          if isrow(dx)
            dx  =  dx';
          end
          if numel(dx) == 1
            dx = dx * ones(size(y));
          end
        otherwise
          disp(['Unknown option ''' num2str(varargin{jj}) ''''])
      end
      jj = jj+1;
    end
  end
  
  % Prepare the vector of weights
  if noerrors
    wt = ones(size(x));
  else
    wt = 1./(dy.^2);
  end
  
  % Compose the ingredients for the textbook LSQ solution for linear model
  Sw  = sum(wt);
  Sx  = sum(x .* wt);
  Sx2 = sum(x.^2 .* wt);
  Sy  = sum(y .* wt);
  sxy = sum(x .* y .* wt);
  
  Delta = Sw*Sx2 - Sx.^2;
  
  %% Calculate the linear fit coefficients using the textbook LSQ solution for linear model
  m = (Sw * sxy - Sx * Sy) / Delta;
  b = (Sx2 * Sy - Sx * sxy) / Delta;
  
  %% Calculate the number of degrees of freedom
  % -2 is for 2 degrees of freedom removed by fit, m and b
  dof = N - 2;
  
  if noerrors
    % The user did not provide dy errors
    % assign the reduced chi^2 a value of 1
    chi2 = 1;
    
    % estimate 'a posteriori' the data uncertainty
    dy_post = sqrt(sum(1 / dof .* (y - (m.*x + b)).^2)) .* ones(size(y));
    
    % call the routine using these estimated uncertainties
    out = regressione_lineare(x, y, dy_post);
    m  = out.m;
    b  = out.b;
    dm = out.dm;
    db = out.db;
    C  = out.C;
    
  else
    % The user did provide some dy errors
    dm = sqrt(Sw / Delta);
    db = sqrt(Sx2 / Delta);
    
    % estimate the reduced chi^2
    chi2 = sum((y - (m*x + b)).^2 .* wt) / dof;
    
    % calculate the linear fit coefficients covariance matrix
    C = (1 / Delta) .* [Sx2 -Sx; -Sx Sw];

    % estimate the contribution of the uncertainties along x, if any
    if ~isempty(dx)
      % propaagate the dx errors via the first parameter estimates
      dy_est = sqrt(dy.^2 + (m .* dx).^2);
      
      % call the routine using these updated uncertainties
      out = regressione_lineare(x, y, dy_est);
      m    = out.m;
      b    = out.b;
      dm   = out.dm;
      db   = out.db;
      C    = out.C;
      chi2 = out.chi2;
      
    end
    
  end
  
  % Prepare the outputs
  out = struct(...
    'm', m, ...
    'b', b, ...
    'dm', dm, ...
    'db', db, ...
    'C',  C, ...
    'chi2', chi2, ...
    'dof',  dof ...
    );
  
end
