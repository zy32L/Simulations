%% P81 奇异函数
%this is one change
clear all;
clc;

% exp(-2t)u(t)
tspan = 0:0.01:2;
y = exp(-2*tspan);
figure;
plot(tspan, y);
ylabel('h(t)=exp(-2t)u(t)');

% 参数定义
tspace = 0:0.0001:2;
t = 0:0.0001:2;
delta_index = [0.25 0.1 0.0025];
DtConvRt = zeros(3,length(tspace)); %便于delta(t)*r(t)
RtConvRt = zeros(3,length(tspace)); %便于r(t)*r(t)

%% delta(t)
figure;
for i = 1:3
    d = zeros(1,length(tspace));
    d(1:delta_index(i)/0.0001) = 1/delta_index(i);
    DtConvRt(i,1:end) = d; %存储delta(t)

    d_interp = @(t) interp1(tspace, d, t);
    
    dydt = @(t, y) -2*y + d_interp(t);
    y0 = 0;
    [t0,y] = ode45(dydt, [0,2], y0);

    plot(t0,y);
    hold on;
end
title('x(t)=\delta(t)');

%% r(t)
figure;
for i = 1:3
    d = zeros(1,length(tspace));
    impulse_num = delta_index(i)/0.0001;

    for j = 1:impulse_num
        d(j) = j/impulse_num*(1/delta_index(i));
    end

    for j = impulse_num+1:impulse_num*2
        d(j) = (2*impulse_num-j)/impulse_num*(1/delta_index(i));
    end
    
    Conv = conv(DtConvRt(i,1:end),d);
    DtConvRt(i,1:end) = Conv(1:length(t)); %delta(t)*r(t)
    Conv = conv(d,d);
    RtConvRt(i,1:end) = Conv(1:length(t)); %r(t)*r(t)

    d_interp = @(t) interp1(tspace, d, t);

    dydt = @(t, y) -2*y + d_interp(t);
    y0 = 0;
    [t0,y] = ode45(dydt, [0,2], y0);

    plot(t0,y);
    hold on;
end
title('x(t)=r(t)');

%% delta(t)*r(t)
figure;
for i = 1:3
    d = DtConvRt(i,1:end);

    d_interp = @(t) interp1(tspace, d, t);

    dydt = @(t, y) -2*y + d_interp(t);
    y0 = 0;
    [t0,y] = ode45(dydt, [0,2], y0);

    plot(t0,y);
    hold on;
end
title('x(t)=\delta(t)*r(t)');

%% r(t)*r(t)
figure;
for i=1:3
    d = RtConvRt(i,1:end);

    d_interp = @(t) interp1(tspace, d, t);

    dydt = @(t, y) -2*y + d_interp(t);
    y0 = 0;
    [t0,y] = ode45(dydt, [0,2], y0);

    plot(t0,y);
    hold on;
end
title('x(t)=\delta(t)*\delta(t)');


