function deltaVec = timeToConv(xvalue,tvalue)

global delta z1Star

z1delta = 0;
timeToDeltaIdx = 1;

% Finding time of convergence for IC: 
    for i=2:length(xvalue(:,1))
        if (((abs(z1Star - xvalue(i,1)) <= delta) && (abs(z1Star - xvalue(i-1,1)) > delta)))
            timeToDeltaIdx = i;
            z1delta = xvalue(i,1);
        end
    end
    deltaVec(1) = z1delta;
    deltaVec(2) = xvalue(timeToDeltaIdx,2);
    deltaVec(3) = tvalue(timeToDeltaIdx,1); 
end