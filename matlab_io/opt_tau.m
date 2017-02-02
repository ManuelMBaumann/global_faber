function tau = opt_tau(e,w,W)

r  = sqrt( w*W*(1.0+e^2) );
th = atan( -sqrt( (e^2*(W+w)^2+(W-w)^2)/(4.0*w*W) ) );
tau = r*cos(th) + 1i*(r*sin(th));

end
