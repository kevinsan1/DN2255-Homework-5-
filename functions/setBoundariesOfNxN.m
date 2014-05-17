function b = setBoundariesOf(myMatrix,sqr)
[sx,sy] = size(myMatrix);
if sqr == 1
    b = zeros(sx+2,sy+2);
    b(2:sx+1,2:sy+1) = myMatrix;
    b(1,:) = b(2,:);
    b(:,1) = b(:,2);
    b(end,:) = b(end-1,:);
    b(:,end)=b(:,end-1);
else
    b = zeros(1,sy+2);
    b(2:sy+1) = myMatrix;
    b(1) = b(2);
    b(end) = b(end-1);
end

