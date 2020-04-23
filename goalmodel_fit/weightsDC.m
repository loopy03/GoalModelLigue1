function [ weights ] = weightsDC( dates , varargin )
% WEIGHTSDC Compute Dixon-Coles weights
%
% Inputs :
%   dates  n by 1 datetime vector of dates
%   Optional (Name,Value) :
%     currentDate  (default : last date in dates) The date which to count
%       backwards from (type datetime)
%     xi  dumping factor (between 0.001 and 0.003, default : 0)
%
% Output :
%   weights  n by 1 double that determine the influence of each match 
%     accorind to dates
%
% Defaults for optional parameters
xi = 0 ;
currentDate = max(dates) ;

% Parse optional input parameters
v = 1;
while v < numel(varargin)
  switch lower(varargin{v})
  case 'currentdate'
    assert(v+1<=numel(varargin));
    v = v+1;
    currentDate = varargin{v};
  case 'xi'
    assert(v+1<=numel(varargin));
    v = v+1;
    xi = varargin{v};  
  otherwise
    error('Unsupported optional input parameter: %s',varargin{v});
  end
  v = v+1;
end

% Compute de weights vector
datediffs = days((dates-currentDate)*-1) ;
weights = exp(-xi*datediffs) ;

% Fixed futur dates to 0
t = find(weights>1,1) ;
weights(t:end) = 0 ;

end

