%% CROSSENTROPY 
% crossentropy comptues the cross-entropy between a true (or posteriori)
% distribution p and an approximation (or prior) distribution q
% 
function E = crossentropy(varargin)

if nargin == 1
    % entropy
    p = ksdensity(varargin{1});
    p = p ./ sum(p);
    p(p==0) = [];
    p(p==0) = [];
    
    E = -p * log2(p)';
    
elseif nargin == 2
    % cross-entropy
    pts = linspace(min([varargin{1}(:);varargin{2}(:)]),max([varargin{1}(:);varargin{2}(:)]),100);
    p = ksdensity(varargin{1},pts);
    q = ksdensity(varargin{2},pts);
    p = p ./ sum(p);
    q = q ./ sum(q);
    id = find(p==0 | q==0);
    p(id) = [];
    q(id) = [];
    
    E = -p * log2(q)';
    
else
    error('not right number of input')
end


end