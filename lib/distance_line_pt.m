%Compute the distance between a line and a point
%
%INPUTS
%2 x 2 matrix of points of the line
%2 x 1 point
%
%OUTPUT
%distance line-point
function dist=distance_line_pt(line,pt)

dist=((line(1,2)-line(1,1))*(line(2,1)-pt(2))-(line(1,1)-pt(1))*(line(2,2)-line(2,1)))...
    /(sqrt((line(1,2)-line(1,1))^2+(line(2,2)-line(2,1))^2));