clear all
close all
clf

K = mmread('K.mtx');
C = mmread('C.mtx');
M = mmread('M.mtx');

n = size(K,1);

damping = 0.5;
f   = linspace(1,10,3);
om  = 2*pi*f;
om  = (1-damping*1i)*om;
tau = opt_tau(damping, real(om(1)), real(om(end)));
eta = om./(om-tau);

I = speye(n);
II = speye(2*n);
O = zeros(n,n);
A = [1i*C K;I O];
B = [M O;O I];

spy(A)
% spy(B)

color=['r','b','g'];


for i=1:length(f)
   mat = A*inv(A-tau*B) - eta(i)*II;  
   figure(9)
   hold on
   e = eig(full(mat));
   plot(real(e),imag(e),strcat('x',color(i)));
   plot_circle(tau,eta(i),9,color(i));
   axis equal
end
