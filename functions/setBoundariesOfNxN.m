function x = setBoundariesOfNxN(myMatrix,method)
% method = nothing for periodic

if nargin < 2
    method = 'periodic BC';
end
switch method
    case 'periodic BC' % periodic boundary conditions
        [sx,sy] = size(myMatrix);
        x = zeros(sx+2,sy+2);
        x(2:sx+1,2:sy+1) = myMatrix;
        x(1,:) = x(end-1,:);
        x(:,1) = x(:,end-1);
        x(end,:) = x(2,:);
        x(:,end)=x(:,2);
    case 'no flux BC' % no flux boundary conditions
        [sx,sy] = size(myMatrix);
        x = zeros(sx+2,sy+2);
        x(2:sx+1,2:sy+1) = myMatrix;
        x(1,:) = x(2,:);
        x(:,1) = x(:,2);
        x(end,:) = x(end-1,:);
        x(:,end)=x(:,end-1);
    case 'reflective velocity BC'
        [sx,sy] = size(myMatrix);
        x = zeros(sx+2,sy+2);
        x(2:sx+1,2:sy+1) = myMatrix;
        x(1,:) = -x(2,:);
        x(:,1) = -x(:,2);
        x(end,:) = -x(end-1,:);
        x(:,end) = -x(:,end-1);
%         disp('reflective BC');
end


