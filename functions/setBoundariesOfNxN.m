function b = setBoundariesOfNxN(myMatrix)
lengthMatrix = length(myMatrix);
b = zeros(lengthMatrix+2,lengthMatrix+2);
b(2:lengthMatrix+1,2:lengthMatrix+1) = myMatrix;
b(1,:) = b(2,:);
b(:,1) = b(:,2);
b(end,:) = b(end-1,:);
b(:,end)=b(:,end-1);