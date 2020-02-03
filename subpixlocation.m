% subpixlocation(line,c,'parabolic','maximum');
function c = subpixlocation(line,c,ignore1,ignore2)

try

d = double(line(c+(-1:1)));
m = d(1)-2*d(2)+d(3);
if m~=0
   c = c + (d(1)-d(3))/(2*m);
end

catch
    c=NaN;
end


