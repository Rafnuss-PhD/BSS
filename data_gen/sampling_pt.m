function pt = sampling_pt(grid,field,method)
% SAMPLING_PT simulate the measurement of high resolution point of the a
% field. Two method are implemented (borehole or random)
% INPUT:
%       - grid:     grid of the matrix to generate (see Import_dat for more info)
%       - field:   true g field
%       - method:   choice of method : 1. borehole | 2. random generated
%       - plotit:   1 or 0 to disply a plot or not
%
% OUTPUT:
%       - pt:        points measurement (pt.x, pt.y, pt.d)
%
% Author: Raphael Nussbaumer
% date : January 2014
% need to do : add assert() for input, flexible number of input var

if isempty(field)
    warning('field is empty')
    pt=[];
    return
end

switch method
    case 1 % select all point in k boreholes equaly spaced
        k=3; % number of borehole;
        %  positions of conditioning data
        
        bor_pos_x = ceil(grid.nx/k/2:grid.nx/k:grid.nx);
        bor_pos_y = 1:grid.ny;
        [x_x, x_y] = meshgrid(grid.x(bor_pos_x),grid.y(bor_pos_y));
        
        % Create the input data
        x = field(bor_pos_y,bor_pos_x);
        pt.d=x(:); pt.x=x_x(:); pt.y=x_y(:); 

    case 2 % select k random point on the mesh
        k= min([500 round(0.01*numel(field(:)))]); % number of input data
        rng(123456)
        [pt.d,idx] = datasample(field(:),k,'Replace',false);
        [pt_j,pt_i] = ind2sub(size(field),idx);
        pt.x=grid.x(pt_i);
        pt.y=grid.y(pt_j);

end
pt.n=numel(pt.d);
end