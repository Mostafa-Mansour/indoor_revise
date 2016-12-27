function [ y] = getRange_( x1,x2,k,wallNum,D )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch wallNum
    case 1
        y=x2/cosd(k-180);
    case 2
        y=x1/cosd(270-k);
    case 3
        y=(D(2)-x2)/cosd(k);
    case 4
        y=(D(1)-x1)/cosd(90-k);
    otherwise
             disp('Pls, check wall Number')
             y=nan;
             
end

end

