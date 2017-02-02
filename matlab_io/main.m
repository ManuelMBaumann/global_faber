clear all
close all
clf

K = mmread('K.mtx');
C = mmread('C.mtx');
M = mmread('M.mtx');


damping = 0.01;
f  = linspace(1,10,5);
om = 2*pi*f;
om = (1-damping*1i)*om;

tau = opt_tau(damping, real(om[1]), real(om[end]));

mat = K + 1i*om[1]*C - om[1]^2*M;
spy(mat); % Have fun...
