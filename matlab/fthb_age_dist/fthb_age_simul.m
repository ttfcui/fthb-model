function fthb_age_simul(p_H, theta, mu, sigma)

    dt = 1/1000;
    T = 40;
    w = zeros(1+T/dt, 1);
    dB_t = sigma*sqrt(dt)*randn(T/dt, 1);

    for i = 1:T/dt
        w(i+1) = w(i) + mu*dt + dB_t(i);
    end

    figure; hold on
    plot(20+ [0:dt:T], w)
    plot(20+ [0:dt:T], ones(1, T/dt+1)*log(p_H*theta), '--k')
    print('-dpdf', 'Brown')
end
