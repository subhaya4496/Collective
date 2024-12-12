% wall interaction
function [FxWall, FyWall, TzWall] = wall_E(d1, d2, FrLS, N, sig)
FxWall = zeros(1,N);
FyWall = zeros(1,N);
TzWall = zeros(1,N);
for cell = 1: N
%     if (abs(d1(cell))< sig)
%         FxWall(cell)= -sign(d1(cell)).* FrLS *(abs(d1(cell)) - sig/2); 
%     end
    if (abs(d2(cell))< sig)
        FyWall(cell)= -sign(d2(cell))* FrLS *(abs(d2(cell)) - sig/2);
    end
end