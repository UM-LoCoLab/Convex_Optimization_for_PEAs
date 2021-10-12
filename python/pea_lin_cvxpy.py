import cvxpy as cp
import numpy as np
from scipy import integrate
from scipy.fft import fft
import math
import scipy.io as sio
import gurobipy
import mosek
import matplotlib.pyplot as plt

def get_a_num(t, dtau, ddtau):
  return integrate.simpson((1/(K_m*K_m))*(-r*(I_m*ddtau+b_m*dtau))*(-r*(I_m*ddtau+b_m*dtau)) + b_m*r*r*dtau*dtau, t)

def get_b_num(t, dq, ddq, tau, dtau, ddtau):
  return integrate.simpson((2/(K_m*K_m))*(-r*(I_m*ddtau+b_m*dtau))*(I_m*r*ddq+b_m*r*dq-(1/r)*tau) - 2*b_m*r*r*dq*dtau, t)

def get_c_num(t, dq, ddq, tau):
  return integrate.simpson((1/(K_m*K_m))*(I_m*r*ddq+b_m*r*dq-(1/r)*tau)*(I_m*r*ddq+b_m*r*dq-(1/r)*tau) + b_m*r*r*dq*dq - tau*dq, t)

def get_coef_SEA(t, dq, ddq, tau, dtau, ddtau):
  """
  Returns the a, b, and c coefficients for SEA optimization using the definition from (8) in ICORR19 paper
  """
  gam1 = -(I_m*ddtau*r + b_m*dtau*r)
  gam2 = I_m*ddq*r + b_m*dq*r - tau/r
  a = np.trapz( gam1**2/K_m**2 + b_m*r**2*dtau**2, t)
  b = np.trapz( 2*gam1*gam2/K_m**2 - 2*b_m*r**2*dq*dtau, t)
  c = np.trapz( gam2**2/K_m**2 + b_m*dq**2*r**2 - dq*tau, t)
  return a, b, c

def get_Q_num(t, q, dtau, ddtau):
  Q_a11 = integrate.simpson((-(I_m*r*ddtau+b_m*r*dtau))*(-(I_m*r*ddtau+b_m*r*dtau)), t)
  Q_a12 = integrate.simpson((-(I_m*r*ddtau+b_m*r*dtau))*((1/r)*q), t)
  Q_a21 = Q_a12
  Q_a22 = integrate.simpson(((1/r)*q)*((1/r)*q), t)
  Q_c11 = integrate.simpson((-r*dtau)*(-r*dtau), t)
  Q_c12 = 0
  Q_c21 = 0
  Q_c22 = 0
  return (1/(K_m*K_m))*np.array([[Q_a11, Q_a12],[Q_a21, Q_a22]]) + b_m*np.array([[Q_c11, Q_c12],[Q_c21, Q_c22]])

def get_f_num(t, beta, q, dq, ddq, tau, dtau, ddtau):
  f_a_1 = integrate.simpson((I_m*r*ddq+b_m*r*dq-(1/r)*tau+beta*(I_m*r*ddq+b_m*r*dq))*(-(I_m*r*ddtau+b_m*r*dtau)), t)
  f_a_2 = integrate.simpson((I_m*r*ddq+b_m*r*dq-(1/r)*tau+beta*(I_m*r*ddq+b_m*r*dq))*((1/r)*q), t)
  f_c_1 = integrate.simpson((r*dq+beta*(r*dq))*(-r*dtau), t)
  f_c_2 = 0
  return (2/(K_m*K_m))*np.array([[f_a_1], [f_a_2]]) + 2*b_m*np.array([[f_c_1], [f_c_2]])

def get_g_num(t, beta, dq, ddq, tau):
  return integrate.simpson((1/(K_m*K_m))*(I_m*r*ddq+b_m*r*dq-(1/r)*tau+beta*(I_m*r*ddq+b_m*r*dq))*(I_m*r*ddq+b_m*r*dq-(1/r)*tau+beta*(I_m*r*ddq+b_m*r*dq))+b_m*(r*dq+beta*(r*dq))*(r*dq+beta*(r*dq))-tau*dq, t)

def FitTrigPoly(output, input, order):
    """
    Fit trigonometric polynomial to the values of `output` as a function of `input`. Order defines the order for the
     sin/cos.

    The function assumes that `input` is periodic when `output` goes from 0-1.
    """
    input = input/input[-1]       # Normalizing data points from (0-1)
    #------------------------------------------- Generate regressor matrix A
    n = int(1 + order*2)
    m = input.size
    A = np.empty([m, n])     # 1 constant term plus sin and cos for each order
    
    A[:, 0] = np.ones(m)
    for i in range(1, order+1):
        A[:, i*2-1] = np.cos(2*np.pi*i*input)
        A[:, i*2]   = np.sin(2*np.pi*i*input)
    
    #------------------------------------------- Fit model to data using cvxpy
    x = cp.Variable(n)
    constraints = []
    # Form objective.
    obj = cp.Minimize( cp.norm(output - A@x, p = 2) )
    # Form and solve problem.
    prob = cp.Problem(obj, constraints)
    prob.solve(solver=cp.MOSEK)  # Returns optimal value  
    if (prob.status == 'optimal'):
        #-------------------------------------- Generate analytical function
        def fun(phase):
            """
            Analytical expression for the fourier series (input is normalized from 0-1)
            """
            acc = x.value[0]
            for i in range(1, order+1):
                cosIdx = i*2-1
                sinIdx = i*2
                acc = acc + x.value[cosIdx]*np.cos(2*np.pi*i*phase)
                acc = acc + x.value[sinIdx]*np.sin(2*np.pi*i*phase)
            return acc
        return fun(input)
    else:
        raise NameError('Optimization program was not optimal')

def get_harmonics(y, pas):
  size=np.shape(y)
  N=size[0]
  Y=fft(y)
  size=np.shape(Y)
  
  w = np.zeros(N)
  for i in range(0,N):
    w[i]=(2*math.pi/(N*pas))*i
  
  w0 = w[1]
  nc2 = math.trunc(N/2)
  ak = np.zeros(nc2)
  bk = np.zeros(nc2)
  for l in range(2, nc2+2):
    ak[l-2]=np.real(Y[l-1]+Y[N-l+1])/N
    bk[l-2]=-np.imag(Y[l-1]-Y[N-l+1])/N

  a0=2*np.real(Y[0])/N;
  return w0, a0, ak, bk

# Global Variables:---------------GLOBAL VARIABLES CHECKED
I_m = 1.2e-4      # Inertia of the motor [Kg*m^2]
b_m = 0.16e-3     # Drag torque [Nm / (rad/s)]
r = 50            # reduction ratio
K_t = 0.14        # Torque constant motor [Nm/A]
R_m = 0.186       # Terminal resistance [ohms]
K_m = K_t/math.sqrt(R_m)    # Motor constant
V_m = 36          # Voltage power supply [Volts]
dq_max = V_m/K_t  # Max motor velocity [rad/sec]
tau_m_max = 28.7*K_t        # Peak torque [Nm]
tau_m_RMS = 7.7*K_t         # Max. continuous current [Nm]
delta_s_max = 90*math.pi/180 # Max. spring elongation [rad]

mc_no = 1000
error = 0.2
cadence='normal'
joint='ankle'

mat_fname = 'InputData/level_' + cadence + '_' + joint + '.mat'
mat_contents = sio.loadmat(mat_fname)
t_data = mat_contents['t'][0]
body_mass= mat_contents['m'][0]
pos_dataRaw = mat_contents['position'][0]*(math.pi/180)
trq_dataRaw = mat_contents['torque'][0]*body_mass
omega0_q, a0_q, ak_q, bk_q = get_harmonics(pos_dataRaw, t_data[1]-t_data[0])
omega0_tau, a0_tau, ak_tau, bk_tau = get_harmonics(trq_dataRaw, t_data[1]-t_data[0])
data_size = np.shape(pos_dataRaw)
data_len = data_size[0]
t = np.linspace(0, 2*math.pi/omega0_q, num=data_len)
pos_data = np.zeros(data_len)+a0_q/2
trq_data = np.zeros(data_len)+a0_tau/2
vel_data = np.zeros(data_len)
dtrq_data = np.zeros(data_len)
acc_data = np.zeros(data_len)
ddtrq_data = np.zeros(data_len)

for idx in range(0, 6):
  pos_data = pos_data + ak_q[idx]*np.cos((idx+1)*omega0_q*t) + bk_q[idx]*np.sin((idx+1)*omega0_q*t)
  trq_data = trq_data + ak_tau[idx]*np.cos((idx+1)*omega0_tau*t) + bk_tau[idx]*np.sin((idx+1)*omega0_tau*t)
  vel_data = vel_data - ((idx+1)*omega0_q)*ak_q[idx]*np.sin((idx+1)*omega0_q*t) + ((idx+1)*omega0_q)*bk_q[idx]*np.cos((idx+1)*omega0_q*t)
  dtrq_data = dtrq_data - ((idx+1)*omega0_tau)*ak_tau[idx]*np.sin((idx+1)*omega0_tau*t) + ((idx+1)*omega0_tau)*bk_tau[idx]*np.cos((idx+1)*omega0_tau*t)
  acc_data = acc_data - ((idx+1)*omega0_q)**2*ak_q[idx]*np.cos((idx+1)*omega0_q*t) - ((idx+1)*omega0_q)**2*bk_q[idx]*np.sin((idx+1)*omega0_q*t)
  ddtrq_data = ddtrq_data - ((idx+1)*omega0_tau)**2*ak_tau[idx]*np.cos((idx+1)*omega0_tau*t) - ((idx+1)*omega0_tau)**2*bk_tau[idx]*np.sin((idx+1)*omega0_tau*t)
  

vel_dataRaw = np.gradient(pos_dataRaw, t)
acc_dataRaw = np.gradient(vel_dataRaw, t)
dtrq_dataRaw = np.gradient(trq_dataRaw, t)
ddtrq_dataRaw = np.gradient(dtrq_dataRaw, t)

# Plotting.
fig, axs = plt.subplots(3, 1)
axs[0].plot(t, pos_dataRaw, label = 'Raw')
axs[0].plot(t, pos_data, '--o', label = 'Fourier')
axs[0].legend()
axs[0].set_xlabel('Time [s]')
axs[0].set_ylabel('Knee Angle [rad] (+) Flexion')
axs[1].plot(t, trq_dataRaw, label = 'Raw - (+) Extension')
axs[1].plot(t, -trq_data, '--o', label = 'Fourier - (+) Flexion')
axs[1].legend()
axs[1].set_xlabel('Time [s]')
axs[1].set_ylabel('Torque [N-m] (+) Extension')
axs[2].plot(t, -vel_data*trq_data)
axs[2].set_xlabel('Time [s]')
axs[2].set_ylabel('Power [W] - (+) Generation')
axs[2].grid()
axs[2].text(0.2, -20, f"Energy per cycle (Provided from the motor to the load): {np.trapz(-vel_data*trq_data, t):3.2f} [J]")
# plt.show()

fig, axs = plt.subplots(3, 1)
axs[0].plot(t, pos_dataRaw, label = 'Raw')
axs[0].plot(t, pos_data, '--', label = 'Fourier')
axs[0].legend()
axs[0].set_xlabel('Time [s]')
axs[0].set_ylabel('Angle [rad] (+) Flexion')
axs[0].grid()
axs[1].plot(t, vel_dataRaw, label = 'Raw')
axs[1].plot(t, vel_data, '--', label = 'Fourier')
axs[1].legend()
axs[1].set_xlabel('Time [s]')
axs[1].set_ylabel('Angular Velocity [rad/s] (+) Flexion')
axs[1].grid()
axs[2].plot(t, acc_dataRaw, label = 'Raw')
axs[2].plot(t, acc_data, '--', label = 'Fourier')
axs[2].legend()
axs[2].set_xlabel('Time [s]')
axs[2].set_ylabel('Angular Acceleration [rad/s^2] (+) Flexion')
axs[2].grid()

fig, axs = plt.subplots(3, 1)
axs[0].plot(t, trq_dataRaw, label = 'Raw')
axs[0].plot(t, trq_data, '--', label = 'Fourier')
axs[0].legend()
axs[0].set_xlabel('Time [s]')
axs[0].set_ylabel('Torque [N-m] (+) Extension')
axs[0].grid()
axs[1].plot(t, dtrq_dataRaw, label = 'Raw')
axs[1].plot(t, dtrq_data, '--', label = 'Fourier')
axs[1].legend()
axs[1].set_xlabel('Time [s]')
axs[1].set_ylabel('Power [N-m/s] (+) Extension')
axs[1].grid()
axs[2].plot(t, ddtrq_dataRaw, label = 'Raw')
axs[2].plot(t, ddtrq_data, '--', label = 'Fourier')
axs[2].legend()
axs[2].set_xlabel('Time [s]')
axs[2].set_ylabel('ddTrq/dt^2 [N-m/s^2] (+) Extension')
axs[2].grid()
# plt.show()


gamma_1 = -r*(I_m*ddtrq_data+b_m*dtrq_data)
gamma_2 = I_m*r*acc_data+b_m*r*vel_data-(1/r)*trq_data

tau_m_0 = (I_m*r*acc_data + b_m*r*vel_data - (1/r)*trq_data)
dq_m_0 = r*vel_data
Em_0 = integrate.simpson((1/(K_m*K_m))*tau_m_0**2 + tau_m_0*dq_m_0, t)

a_num = get_a_num(t, dtrq_data, ddtrq_data)
b_num = get_b_num(t, vel_data, acc_data, trq_data, dtrq_data, ddtrq_data)
c_num = get_c_num(t, vel_data, acc_data, trq_data)
a_test, b_test, c_test = get_coef_SEA(t, vel_data, acc_data, trq_data, dtrq_data, ddtrq_data)
feas_judge = 0


# Linear PEA Optimization
# Create scalar optimization variable - k_p:
k = cp.Variable()

# Create constraints.
constraints = [k >= 0]
a_RMS = np.zeros((2,2))
b_RMS = np.zeros((1,2))
c_RMS = 0
beta = 0

for i in range(len(pos_data)):
  a_T=np.array([[-(I_m*r*ddtrq_data[i]+b_m*r*dtrq_data[i]), (1/r)*pos_data[i]]]) # [:, np.newaxis]
  b=I_m*r*acc_data[i] + b_m*r*vel_data[i] - (1/r)*trq_data[i] + beta*(I_m*r*acc_data[i] + b_m*r*vel_data[i])
  c_T=np.array([[-r*dtrq_data[i], 0]])
  d=r*vel_data[i] + beta*(r*vel_data[i])
  constraints.append(a_T @ cp.hstack([0,k]) + b <= tau_m_max)
  constraints.append(-a_T @ cp.hstack([0,k]) - b <= tau_m_max)
  constraints.append((a_T+(K_m*K_m)*c_T) @ cp.hstack([0,k]) + (b+K_m*K_m*d-K_t*V_m/R_m) <= 0)
  constraints.append((-a_T+(K_m*K_m)*c_T) @ cp.hstack([0,k]) + (-b+K_m*K_m*d-K_t*V_m/R_m) <= 0)
  constraints.append((a_T-(K_m*K_m)*c_T) @ cp.hstack([0,k]) + (b-K_m*K_m*d-K_t*V_m/R_m) <= 0)
  constraints.append((-a_T-(K_m*K_m)*c_T) @ cp.hstack([0,k]) - (b+K_m*K_m*d+K_t*V_m/R_m) <= 0)
  a_RMS += np.transpose(a_T) @ a_T
  b_RMS += (2*b) * a_T
  c_RMS += b*b

n = len(trq_data)
a_RMS = (1/n)*a_RMS
b_RMS = (1/n)*b_RMS
c_RMS = (1/n)*c_RMS

constraints.append((cp.quad_form(cp.hstack([0,k]), a_RMS) + b_RMS @ cp.hstack([0,k]) + c_RMS) <= tau_m_RMS*tau_m_RMS)

Q = get_Q_num(t, pos_data, dtrq_data, ddtrq_data)
f = get_f_num(t, beta, pos_data, vel_data, acc_data, trq_data, dtrq_data, ddtrq_data)
g = get_g_num(t, beta, vel_data, acc_data, trq_data)
# Form objective.
obj = cp.Minimize(cp.quad_form(cp.hstack([0,k]), Q) + f.T @ cp.hstack([0,k]) + g)

# Form and solve problem.
prob = cp.Problem(obj, constraints)
print(f"2nd: PEA Optimization")
print(f"Problem is DCP: {prob.is_dcp()}")
prob.solve(solver=cp.MOSEK)
print("Status:", prob.status)
if prob.status == 'infeasible':
  feas_judge = 1
else:
  tau_m_opt = (I_m*r*acc_data+b_m*r*vel_data)*0*k.value - (I_m*r*ddtrq_data+b_m*r*dtrq_data)*0 + (1/r)*k.value*pos_data + (I_m*r*acc_data + b_m*r*vel_data - (1/r)*trq_data)
  dq_m_opt = r*0*k.value*vel_data - r*0*dtrq_data + r*vel_data
  Em_true = integrate.simpson((1/(K_m*K_m))*tau_m_opt**2 + tau_m_opt*dq_m_opt, t)
  Em_per = (Em_true-Em_0)*100/Em_0
  print(f"k_p: {k.value:.6f}; Em_obj: {prob.value:.6f}; Em_true: {Em_true:.6f} ({Em_per:.1f}%)")

  pos_data_opt_rms = np.sqrt(np.mean(pos_data_opt**2))
  vel_data_opt_rms = np.sqrt(np.mean(vel_data_opt**2))
  acc_data_opt_rms = np.sqrt(np.mean(acc_data_opt**2))
  trq_data_opt_rms = np.sqrt(np.mean(trq_data_opt**2))

  Em_new = np.zeros(mc_no)
  for j in range(0, mc_no):
    # # Noise w.r.t. RMS value:
    pos_data_new = pos_data_opt + pos_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    vel_data_new = vel_data_opt + vel_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    acc_data_new = acc_data_opt + acc_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    trq_data_new = trq_data_opt + trq_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    
    tau_m_new = (1/r)*k.value*pos_data_new + (I_m*r*acc_data_new + b_m*r*vel_data_new - (1/r)*trq_data_new)
    dq_m_new = r*vel_data_new
    Em_new[j] = integrate.simpson((1/(K_m*K_m))*tau_m_new**2 + tau_m_new*dq_m_new, t)


  Em_trial=np.mean(Em_new)
  Em_std = np.std(Em_new)

  Em_act_per = (Em_trial-Em_0)*100/Em_0
  Em_std_act = (Em_trial+2*Em_std-Em_0)*100/Em_0 - Em_act_per
  print(f"Em_trial: {Em_act_per:.1f}%; stdev: {Em_std_act:.1f}%")