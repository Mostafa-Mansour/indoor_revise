function [ y] = getMesurement( x1,x2,k,wallNum,D )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch wallNum
    case 1
        y=x2/cosd(k);
    case 2
        y=x1/sind(k);
    case 3
        y=(D(2)-x2)/cosd(k);
    case 4
        y=(D(1)-x1)/sind(k);
    otherwise
             disp('Pls, check wall Number')
             y=nan;
             
end

end

