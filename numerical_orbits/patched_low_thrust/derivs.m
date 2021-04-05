function dydx = derivs(x,y,nu,itan)
%--------------------------------------------------------------------------
%...spiral orbit equations of motion
if itan == 0  %...circumferential thrust
    dydx(1)=y(2)*y(3)*y(3)-1.0/(y(2)*y(2));
    dydx(2)=y(1);
    dydx(3)=-2.0*y(1)*y(3)/y(2)+nu/y(2);
    dydx(4)=y(3);
end
if itan == -1  %...tangential thrust
    if y(1) > 0
        phi=atan(y(2)*y(3)/y(1));
    end
    if y(1) < 0
        phi=atan(y(2)*y(3)/y(1))+pi;
    end
    if y(1) == 0
        phi=pi/2.0;
    end
    dydx(1)=y(2)*y(3)*y(3)-1.0/(y(2)*y(2))+nu*cos(phi);
    dydx(2)=y(1);
    dydx(3)=-2.0*y(1)*y(3)/y(2)+nu*sin(phi)/y(2);
    dydx(4)=y(3);
end
if itan > 0  %...vectored thrust
    phi=itan*pi/180;
    dydx(1)=y(2)*y(3)*y(3)-1.0/(y(2)*y(2))+nu*cos(phi);
    dydx(2)=y(1);
    dydx(3)=-2.0*y(1)*y(3)/y(2)+nu*sin(phi)/y(2);
    dydx(4)=y(3);
return
end
