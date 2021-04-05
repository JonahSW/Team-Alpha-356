function yout = rk4(y,dydx,n,x,h,nu,itan)
%--------------------------------------------------------------------------
%...4th order Runge-Kutta numerical integration
%...taken from 'Numerical Recipies in FORTRAN' 2nd Edition, page 706
hh=h*0.5;
h6=h/6.0;
xh=x+hh;
%...first step
for i=1:n
    yt(i)=y(i)+hh*dydx(i);
end
%...second step
dyt=derivs(xh,yt,nu,itan);
for i=1:n
    yt(i)=y(i)+hh*dyt(i);
end
%...third step
dym=derivs(xh,yt,nu,itan);
for i=1:n
    yt(i)=y(i)+h*dym(i);
    dym(i)=dyt(i)+dym(i);
end
%...fourth step
dyt=derivs(x+h,yt,nu,itan);
for i=1:n
    yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0*dym(i));
end
return
end
