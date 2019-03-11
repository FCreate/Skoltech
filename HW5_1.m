n = 1000;
t = (0:n)';
dev = 250; 
temp = ones(dev,1);
figure(1)
subplot(311)
%%Signals 
%piecewise linear continious
lc_signal= [t(1:dev); 2*dev-t(dev+1:2*dev);t(2*dev+1:3*dev)-2*dev; 4*dev-t(3*dev+1:4*dev+1)];
lc_signal = lc_signal/(max(lc_signal));
noise = 0.1*randn(size(t));
lc_signal_noise = noise+lc_signal;
plot(t, lc_signal_noise, '-');
xlabel('x')
ylabel('y')
title('Piecewise linear continious')

subplot(312)
%piecewise linear non-continious
lcnc_signal= [t(1:dev); dev-t(dev+1:2*dev);t(2*dev+1:3*dev)-2*dev; dev-t(3*dev+1:4*dev+1)];
lcnc_signal = lcnc_signal/(max(lcnc_signal));
noise = 0.1*randn(size(t));
lcnc_signal_noise = noise+lcnc_signal;
plot(t, lcnc_signal_noise, '-');
xlabel('x')
ylabel('y')
title('Piecewise linear non-continious')

subplot(313)
%sin function
sin_signal= sin((pi/100)*t);
noise = 0.2*randn(size(t));
sin_signal_noise = noise+sin_signal;
plot(t, sin_signal_noise, '-');
xlabel('x')
ylabel('y')
title('Sin function')

%Let's took gamma = 1 and run recovery function for all signals
gamma = 1;
D1 = sparse(n,n+1);
D1(:,1:n) = speye(n,n); 
D1(:,2:n+1) = D1(:,2:n+1)-speye(n,n);
D2 = zeros(n,n+1);
for v = 1:n-1
   D2(v,v) = 1;
   D2(v, v+1)=-2;
   D2(v, v+2) = 1;
end

figure(2)
%Recovery first signal
%D1 norm 2
subplot(311)
cvx_begin quiet
    variable lc_d1_2(n+1)
    minimize(norm(lc_d1_2-lc_signal_noise)+gamma*norm(D1*lc_d1_2))
cvx_end
disp(norm(lc_d1_2 - lc_signal_noise));
plot(t, lc_d1_2, '-');
title('D1 norm 2');

%D1 norm 1
subplot(312)
cvx_begin quiet
    variable lc_d1_1(n+1)
    minimize(norm(lc_d1_1-lc_signal_noise)+gamma*norm(D1*lc_d1_1, 1))
cvx_end
plot(t, lc_d1_1, '-');
title('D1 norm 1');

%D2 norm 1
subplot(313)
cvx_begin quiet
    variable lc_d2_1(n+1)
    minimize(norm(lc_d2_1-lc_signal_noise)+gamma*norm(D2*lc_d2_1, 1))
cvx_end
plot(t, lc_d2_1, '-');
title('D2 norm 1');


figure(3)
%Recovery second signal
%D1 norm 2
subplot(311)
cvx_begin quiet
    variable lc_d1_2(n+1)
    minimize(norm(lc_d1_2-lcnc_signal_noise)+gamma*norm(D1*lc_d1_2))
cvx_end
plot(t, lc_d1_2, '-');
title('D1 norm 2');

%D1 norm 1
subplot(312)
cvx_begin quiet
    variable lc_d1_1(n+1)
    minimize(norm(lc_d1_1-lcnc_signal_noise)+gamma*norm(D1*lc_d1_1, 1))
cvx_end
plot(t, lc_d1_1, '-');
title('D1 norm 1');

%D2 norm 1
subplot(313)
cvx_begin quiet
    variable lc_d2_1(n+1)
    minimize(norm(lc_d2_1-lcnc_signal_noise)+gamma*norm(D2*lc_d2_1, 1))
cvx_end
plot(t, lc_d2_1, '-');
title('D2 norm 1');


figure(4)
%Recovery third signal
%D1 norm 2
subplot(311)
cvx_begin quiet
    variable lc_d1_2(n+1)
    minimize(norm(lc_d1_2-sin_signal_noise)+gamma*norm(D1*lc_d1_2))
cvx_end
plot(t, lc_d1_2, '-');
title('D1 norm 2');

%D1 norm 1
subplot(312)
cvx_begin quiet
    variable lc_d1_1(n+1)
    minimize(norm(lc_d1_1-sin_signal_noise)+gamma*norm(D1*lc_d1_1, 1))
cvx_end
plot(t, lc_d1_1, '-');
title('D1 norm 1');

%D2 norm 1
subplot(313)
cvx_begin quiet
    variable lc_d2_1(n+1)
    minimize(norm(lc_d2_1-sin_signal_noise)+gamma*norm(D2*lc_d2_1, 1))
cvx_end
plot(t, lc_d2_1, '-');
title('D2 norm 1');

%Gamma dependency
figure(5)
gamma_var = linspace(0.01, 10, 100);
d12_norms = zeros(1,100);
yx_norms = zeros(1,100);
for i = 1:length(gamma_var)
    cvx_begin quiet
        variable lc_d1_2(n+1)
        minimize(norm(lc_d1_2-lc_signal_noise)+gamma_var(i)*norm(D1*lc_d1_2))
    cvx_end
    yx_norms(i) = norm(lc_signal_noise-lc_d1_2);
    d12_norms(i) = norm(D1*lc_d1_2);
end
plot(d12_norms,yx_norms, '*');
ylabel('|y-x|');
xlabel('|D1x|2');
title('Gamma dependency first signal');

figure(6)
gamma_var = linspace(0.01, 10, 100);
d12_norms = zeros(1,100);
yx_norms = zeros(1,100);
for i = 1:length(gamma_var)
    cvx_begin quiet
        variable lc_d1_2(n+1)
        minimize(norm(lc_d1_2-lcnc_signal_noise)+gamma_var(i)*norm(D1*lc_d1_2))
    cvx_end
    yx_norms(i) = norm(lcnc_signal_noise-lc_d1_2);
    d12_norms(i) = norm(D1*lc_d1_2);
end
plot(d12_norms,yx_norms, '*');
ylabel('|y-x|');
xlabel('|D1x|2');
title('Gamma dependency second signal');

figure(7)
gamma_var = linspace(0.01, 10, 100);
d12_norms = zeros(1,100);
yx_norms = zeros(1,100);
for i = 1:length(gamma_var)
    cvx_begin quiet
        variable lc_d1_2(n+1)
        minimize(norm(lc_d1_2-sin_signal_noise)+gamma_var(i)*norm(D1*lc_d1_2))
    cvx_end
    yx_norms(i) = norm(sin_signal_noise-lc_d1_2);
    d12_norms(i) = norm(D1*lc_d1_2);
end
plot(d12_norms,yx_norms, '*');
ylabel('|y-x|');
xlabel('|D1x|2');
title('Gamma dependency sin signal');






