pos_data_rms = rms(Leg2.trajectory.ql);
vel_data_rms = rms(Leg2.trajectory.qld);
acc_data_rms = rms(Leg2.trajectory.qldd);
trq_data_rms = rms(Leg2.trajectory.torque);
dtrq_data_rms = rms(Leg2.trajectory.torqued);
n = length(Leg2.trajectory.ql);
error = 0.2;
e_l = -1;
e_u = 1;
mc_no_set = [10,100,1000,5000,10000,50000,75000];

%% Em (Non-linear)
Em_trial = zeros(length(mc_no_set),1);
Em_std = zeros(length(mc_no_set),1);
for j=1:length(mc_no_set)
    Em_new = zeros(mc_no,1);
    for i=1:mc_no
        pos_data_new = Leg2.trajectory.ql + pos_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        vel_data_new = Leg2.trajectory.qld + vel_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        acc_data_new = Leg2.trajectory.qldd + acc_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        trq_data_new = Leg2.trajectory.torque + trq_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        e_new = Leg2.PEAParam.iM * Leg2.PEAParam.r * acc_data_new + Leg2.PEAParam.bm * Leg2.PEAParam.r * vel_data_new + (1/Leg2.PEAParam.r) * trq_data_new;
        F_new = sparse((1/Leg2.PEAParam.r)*diag(pos_data_new));
        tau_m_new = e_new + F_new * Leg2.nonlinearPEA.compliance;
        dq_m_new = Leg2.PEAParam.r * vel_data_new;
        Em_new(i) = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_new.^2 + tau_m_new.*dq_m_new);
    end
    Em_trial(j) = mean(Em_new);
    Em_std(j) = std(Em_new);
end
tau_m_0 = Leg2.PEAParam.iM * Leg2.PEAParam.r * Leg2.trajectory.qldd + Leg2.PEAParam.bm * Leg2.PEAParam.r * Leg2.trajectory.qld + (1/Leg2.PEAParam.r) * Leg2.trajectory.torque;
dq_m_0 = Leg2.PEAParam.r * Leg2.trajectory.qld;
Em_0 = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_0.^2 + tau_m_0.*dq_m_0);
Em_per = (Leg2.nonlinearPEA.energy_total-Em_0)*100/Em_0;
disp(['NL: Em_0: ' num2str(Em_0) '; Em_nom: ' num2str(Leg2.nonlinearPEA.energy_total) ' (' num2str(Em_per) '%)'])

Em_act_per = (Em_trial(end)-Em_0)*100/Em_0;
Em_std_act = (Em_trial(end)+2*Em_std(end)-Em_0)*100/Em_0 - Em_act_per;
disp(['Em_trial: ' num2str(Em_trial(end)) ' (' num2str(Em_act_per) '%); stdev: ' num2str(Em_std_act) '%'])

figure()
semilogx(mc_no_set, Em_trial,'DisplayName','Monte Carlo')
hold on
semilogx([mc_no_set(1) mc_no_set(end)], [Leg2.nonlinearPEA.energy_total Leg2.nonlinearPEA.energy_total],'k--','DisplayName','Nominal')
semilogx([mc_no_set(1) mc_no_set(end)], [Em_0 Em_0],'k-','DisplayName','$E_{m,0}$')
xlabel('No. of Simulations','Interpreter','latex','FontSize',16)
l1=legend('show');
legend 'boxoff'
set(l1,'Interpreter','latex','FontSize',14)
ylabel('$E_{m}$','Interpreter','latex','FontSize',16)
title('Monte Carlo Sim. Convergence','Interpreter','latex','FontSize',16)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
%% Em (Linear)
Em_trial = zeros(length(mc_no_set),1);
Em_std = zeros(length(mc_no_set),1);
for j=1:length(mc_no_set)
    Em_new = zeros(mc_no,1);
    for i=1:mc_no_set(j)
        pos_data_new = Leg2.trajectory.ql + pos_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        vel_data_new = Leg2.trajectory.qld + vel_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        acc_data_new = Leg2.trajectory.qldd + acc_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        trq_data_new = Leg2.trajectory.torque + trq_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
        e_new = Leg2.PEAParam.iM * Leg2.PEAParam.r * acc_data_new + Leg2.PEAParam.bm * Leg2.PEAParam.r * vel_data_new + (1/Leg2.PEAParam.r) * trq_data_new;
        F_new = sparse((1/Leg2.PEAParam.r)*diag(pos_data_new));
        tau_m_new = e_new + F_new * Leg2.linearPEA.compliance;
        dq_m_new = Leg2.PEAParam.r * vel_data_new;
        Em_new(i) = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_new.^2 + tau_m_new.*dq_m_new);
    end
    Em_trial(j) = mean(Em_new);
    Em_std(j) = std(Em_new);
end
Em_per = (Leg2.linearPEA.energy_total-Em_0)*100/Em_0;
disp(['L: Em_0: ' num2str(Em_0) '; Em_nom: ' num2str(Leg2.linearPEA.energy_total) ' (' num2str(Em_per) '%)'])

Em_act_per = (Em_trial(end)-Em_0)*100/Em_0;
Em_std_act = (Em_trial(end)+2*Em_std(end)-Em_0)*100/Em_0 - Em_act_per;
disp(['Em_trial: ' num2str(Em_trial(end)) ' (' num2str(Em_act_per) '%); stdev: ' num2str(Em_std_act) '%'])

figure()
semilogx(mc_no_set, Em_trial,'DisplayName','Monte Carlo')
hold on
semilogx([mc_no_set(1) mc_no_set(end)], [Em_0 Em_0],'k--','DisplayName','Nominal')
xlabel('No. of Simulations','Interpreter','latex')
l1=legend('show');
set(l1,'Interpreter','latex')
ylabel('$E_{m}$','Interpreter','latex')
title('Monte Carlo Sim. Convergence','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')