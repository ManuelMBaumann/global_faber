function [center,radius] = plot_circle(tau,eta,fignr,color)


NOP = 1000;

center(1) = real((0 - conj(tau))/(tau -conj(tau)) - eta);
center(2) = imag((0 - conj(tau))/(tau -conj(tau)) - eta);
radius = abs((tau - 0)/(tau -conj(tau)));
            
% radius2 = 0.5*sqrt( 1+( real(tau)/imag(tau) )^2 )
% center2 = 0.5+(real(tau)/(2*imag(tau)))*i - eta
 
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);

figure(fignr)
hold on
plot(X,Y,char(strcat(color,'-.')),'LineWidth',2);
plot(center(1),center(2),char(strcat(color,'+')),'MarkerSize',12,'LineWidth',1);
axis equal;

end
