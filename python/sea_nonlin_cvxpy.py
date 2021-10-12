import cvxpy as cp
import numpy as np
from scipy import integrate
from scipy.fft import fft
import math
import scipy.io as sio
import gurobipy
import mosek
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from datetime import datetime
import os

def get_a_num(dq):
  return r*dq

def get_B_num(dtau):
  return -r*np.diag(dtau)

def get_c_num(ddq):
  return r*ddq

def get_D_num(dtau, ddtau):
  data_len = len(dtau)
  D_1 = np.zeros((data_len, data_len))
  for idx in range(0, data_len):
    cidx=(idx+1)%data_len
    D_1[idx][cidx]=1
    cidx=(idx+data_len-1)%data_len
    D_1[idx][cidx]=-1

  D_1 = D_1*(1/(2*dt))
  return -r*np.diag(ddtau)-r*(np.diag(dtau)@D_1)

def get_e_num(a, c, tau):
  return I_m*c+b_m*a-(1/r)*tau

def get_F_num(B, D):
  return I_m*D+b_m*B

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

def get_G_num(B, F):
  return dt*((1/(K_m*K_m))*(F.T @ F) + b_m*(B.T @ B))

def get_h_num(a, B, e, F):
  return dt*((2/(K_m*K_m))*(e.T @ F) + 2*b_m*(a.T @ B))

def get_w_num(a, e, dq, tau):
  return dt*((1/(K_m*K_m))*(e.T @ e) + b_m*(a.T @ a) - tau.T @ dq)

def get_elong_matr(tau, dtau, argidx):
  mintauidx=argidx[0]
  n=len(tau)
  i=np.zeros(len(argidx))
  for idx in range(0, len(argidx)):
    i[idx]=idx

  ordMatr = csr_matrix((np.ones(len(argidx)), (i, argidx)), shape = (n, n)).toarray()
  elo_mat = np.zeros((n,n))
  elo_mat[1][0]=dtau[0];
  elo_mat[1][1]=dtau[1];
  for i in range(2, len(dtau)):
    for j in range(0,n):
      elo_mat[i][j]=elo_mat[i-1][j]
    elo_mat[i][i] = dtau[i]

  for i in range(0,n):
    elo_mat[i][0]=elo_mat[i][0]/2

  for i in range(0,n):
    elo_mat[i][i]=elo_mat[i][i]/2

  del_vector = np.zeros(n)
  for i in range(0,n):
    del_vector[i] = elo_mat[mintauidx][i]

  for i in range(0, len(dtau)):
    for j in range(0, n):
      elo_mat[i][j] = elo_mat[i][j] - del_vector[j]

  elo_mat = dt*elo_mat
  elo_mat = ordMatr@elo_mat
  return elo_mat

def get_L_num(dtau_s):
  data_size = np.shape(dtau_s)
  data_len = data_size[0]
  L = np.zeros((data_len, data_len))
  L[1][0]=dtau_s[0]
  L[1][1]=dtau_s[1]
  for idx in range(2, data_len):
    for cidx in range(0, idx+1):
      if cidx==0 or cidx==idx:
        L[idx][cidx]=dtau_s[cidx]
      else:
        L[idx][cidx]=2*dtau_s[cidx]

  return (dt/2)*L

def get_diff_mat(n, argidx):
  D_M = np.zeros((n-1, n))
  for ridx in range(0, n-1):
    minus_idx = argidx[ridx]
    add_idx = argidx[ridx+1]
    D_M[ridx][minus_idx]=-1
    D_M[ridx][add_idx]=1
  
  return D_M

def get_L_diff_mat(L, n):
  L_star = np.zeros((n-1, n))
  L_pl = np.zeros((n-1, n))
  rcount = 0
  for ridx in range(1, n):
    for cidx in range(0, n):
      L_star[rcount][cidx]=L[ridx][cidx]

    rcount=rcount+1

  rcount = 0
  for ridx in range(0, n-1):
    for cidx in range(0, n):
      L_pl[rcount][cidx]=L[ridx][cidx]

    rcount=rcount+1

  return L_star-L_pl

def get_elo_diff_matr(L_e):
  L_e_f = np.zeros((data_len-2,data_len-1))
  L_e_r = np.zeros((data_len-2,data_len-1))
  for i in range(0, data_len-2):
    for j in range(0, data_len-1):
      L_e_f[i][j]=L_e[i+1][j]
      L_e_r[i][j]=L_e[i][j]
  return -(L_e_f - L_e_r)

def get_spr_int_vec(tau):
  difftau=np.diff(tau)
  int_vec = np.append(difftau,[0])+np.insert(difftau,0,[0])
  int_vec[0] = int_vec[0]/2
  int_vec[-1] = int_vec[-1]/2
  int_vec = int_vec*tau
  return int_vec

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
dt = t[1]-t[0]
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
# axs[1].plot(t, trq_dataRaw, label = 'Raw - (+) Extension')
axs[1].plot(t, -trq_data, '--o', label = 'Fourier (+) Flexion')
axs[1].legend()
axs[1].set_xlabel('Time [s]')
axs[1].set_ylabel('Torque [N-m] (+) Flexion')
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

pos_data_opt = np.zeros(data_len-1)
trq_data_opt = np.zeros(data_len-1)
vel_data_opt = np.zeros(data_len-1)
dtrq_data_opt = np.zeros(data_len-1)
acc_data_opt = np.zeros(data_len-1)
ddtrq_data_opt = np.zeros(data_len-1)

for i in range(0, data_len-1):
  pos_data_opt[i] = pos_data[i]
  trq_data_opt[i] = trq_data[i]
  vel_data_opt[i] = vel_data[i]
  dtrq_data_opt[i] = dtrq_data[i]
  acc_data_opt[i] = acc_data[i]
  ddtrq_data_opt[i] = ddtrq_data[i]

a = get_a_num(vel_data_opt)
B_num = get_B_num(dtrq_data_opt)
c = get_c_num(acc_data_opt)
D_num = get_D_num(dtrq_data_opt, ddtrq_data_opt)
e = get_e_num(a, c, trq_data_opt)
F_num = get_F_num(B_num, D_num)
G_num = get_G_num(B_num, F_num)
h = get_h_num(a, B_num, e, F_num)
w_num = get_w_num(a, e, vel_data_opt, trq_data_opt)
L_num = get_L_num(-dtrq_data_opt)
idxsorted = np.argsort(trq_data_opt)
L_e = get_elong_matr(trq_data_opt, dtrq_data_opt, idxsorted)
diff_mat = get_diff_mat(data_len-1, idxsorted)
L_diff = get_L_diff_mat(L_num, data_len-1)
L_elo_diff = get_elo_diff_matr(L_e)
tau_dtau = (-trq_data)*(-dtrq_data)
L_tau_mat = get_L_num(tau_dtau)
L_tau_vec = get_spr_int_vec(trq_data)
a_num = np.zeros((data_len-1, 1))
c_num = np.zeros((data_len-1, 1))
e_num = np.zeros((data_len-1, 1))
h_num = np.zeros((data_len-1, 1))
for i in range(0, data_len-1):
  a_num[i][0]=a[i]
  c_num[i][0]=c[i]
  e_num[i][0]=e[i]
  h_num[i][0]=h[i]

feas_judge = 0

# Nonlinear SEA Optimization
# Create scalar optimization variable - the compliance:
alpha_s = cp.Variable((data_len-1, 1))

# Create constraints.
constraints = []

# # alpha_s >= 0
constraints.append(diff_mat @ (L_num @ alpha_s) <= np.zeros((data_len-2,1)))

# # E_s=0
constraints.append(L_tau_vec @ cp.vstack((alpha_s,alpha_s[np.newaxis,0])) <= 0)
constraints.append(-L_tau_vec @ cp.vstack((alpha_s,alpha_s[np.newaxis,0])) <= 0)


# ------- Actuator Constraints -------
constraints.append(L_e @ alpha_s <= np.ones((data_len-1, 1))*delta_s_max)
constraints.append(-L_e @ alpha_s <= np.ones((data_len-1, 1))*delta_s_max)
constraints.append(B_num @ alpha_s <= np.ones((data_len-1, 1))*dq_max - a_num)
constraints.append(-B_num @ alpha_s <= np.ones((data_len-1, 1))*dq_max + a_num)

constraints.append(F_num @ alpha_s <= np.ones((data_len-1, 1))*tau_m_max - e_num)
constraints.append(-F_num @ alpha_s <= np.ones((data_len-1, 1))*tau_m_max + e_num)

constraints.append((F_num+(K_t**2/R_m)*B_num) @ alpha_s <= np.ones((data_len-1, 1))*(K_t*V_m/R_m) - e_num - a_num*(K_t**2/R_m))
constraints.append((F_num-(K_t**2/R_m)*B_num) @ alpha_s <= np.ones((data_len-1, 1))*(K_t*V_m/R_m) - e_num + a_num*(K_t**2/R_m))
constraints.append((-F_num+(K_t**2/R_m)*B_num) @ alpha_s <= np.ones((data_len-1, 1))*(K_t*V_m/R_m) + e_num - a_num*(K_t**2/R_m))
constraints.append((-F_num-(K_t**2/R_m)*B_num) @ alpha_s <= np.ones((data_len-1, 1))*(K_t*V_m/R_m) + e_num + a_num*(K_t**2/R_m))

p_1 = np.ones((data_len-1, 1))*dq_max - a_num
p_2 = np.ones((data_len-1, 1))*dq_max + a_num
p_3 = np.ones((data_len-1, 1))*tau_m_max - e_num
p_4 = np.ones((data_len-1, 1))*tau_m_max + e_num
p_5 = np.ones((data_len-1, 1))*(K_t*V_m/R_m) - e_num - a_num*(K_t**2/R_m)
p_6 = np.ones((data_len-1, 1))*(K_t*V_m/R_m) - e_num + a_num*(K_t**2/R_m)
p_7 = np.ones((data_len-1, 1))*(K_t*V_m/R_m) + e_num - a_num*(K_t**2/R_m)
p_8 = np.ones((data_len-1, 1))*(K_t*V_m/R_m) + e_num + a_num*(K_t**2/R_m)

M_num = np.vstack([B_num, -B_num, F_num, -F_num, F_num+(K_t**2/R_m)*B_num, F_num-(K_t**2/R_m)*B_num, -F_num+(K_t**2/R_m)*B_num, -F_num-(K_t**2/R_m)*B_num])
p_num = np.vstack((p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8))

u, s, vh = np.linalg.svd(F_num, full_matrices=True)
# FF = vh.T.dot(np.diag(s)**2).dot(vh)
FF = cp.Parameter((data_len-1, data_len-1), PSD=True)
FF.value = vh.T @ np.diag(s)**2 @ vh
# FF.value = F_num.T @ F_num
# FF = F_num.T @ F_num

eF = e_num.T @ F_num
ee = e_num.T @ e_num
constraints.append(cp.quad_form(alpha_s, FF) + (2*eF)@alpha_s + ee <= (data_len-1)*(tau_m_RMS)**2)

# Form objective.
obj = cp.Minimize(cp.quad_form(alpha_s, G_num) + h_num.T @ alpha_s + w_num)

# Form and solve problem.
prob = cp.Problem(obj, constraints)
print(f"1st: SEA Optimization")
print("Problem is DCP:", prob.is_dcp())
prob.solve(solver=cp.MOSEK)
print("Status:", prob.status)
if prob.status == 'infeasible':
  feas_judge = 1
else:
  
  fig, axs = plt.subplots(1, 1)
  POS = np.zeros((data_len-1, 1))
  TRQ = np.zeros((data_len-1, 1))
  TRQ_s = np.zeros(data_len-1)
  ALP = np.zeros((data_len-1, 1))
  DISP = np.zeros((data_len-1, 1))
  DISPDot = np.zeros((data_len-1, 1))
  dTRQ_s = np.zeros(data_len-1)
  delta_s = np.zeros(data_len-1)
  
  for i in range(0, data_len-1):
    POS[i][0]=pos_data[i]
    TRQ[i][0]=-trq_data[i]
    TRQ_s[i] = -trq_data[i]
    ALP[i][0]=alpha_s.value[i][0]
    dTRQ_s[i]=-dtrq_data[i]
    
  elong = POS-(1/r)*(r*POS-r*(np.diag(-TRQ_s)@ALP))
  DISP[0][0]=POS[0][0]-(1/r)*(r*POS[0][0]-r*(TRQ_s[0]*ALP[0][0]))
  DISPDot[0][0]=vel_data[0]-(1/r)*(r*vel_data[0]-r*((-dtrq_data[0])*ALP[0][0]))
  for i in range(1, data_len-1):
    DISP[i][0]=DISP[i-1][0]+ALP[i-1][0]*(TRQ_s[i]-TRQ_s[i-1])
    DISPDot[i][0]=vel_data[i]-(1/r)*(r*vel_data[i]-r*((-dtrq_data[i])*ALP[i][0]))

  EL = np.zeros(data_len-1)
  for i in range(0, data_len-1):
    EL[i] = DISP[i][0]
  ALP_mask = np.gradient(EL, TRQ_s)
  delta_s[0]=0
  for i in range(1, data_len-1):
    DISPDot_list = np.zeros(i)
    t_set = np.zeros(i)
    for j in range(0, i):  
      DISPDot_list[j]=DISPDot[j][0]
      t_set[j] = t[j]

    delta_s[i]=integrate.simpson(DISPDot_list, t_set)

  ALPHA_vec = np.zeros(data_len)
  for i in range(0, data_len-1):
    ALPHA_vec[i]=alpha_s.value[i][0]

  ALPHA_vec[-1]=alpha_s.value[0][0]
  a_opt = get_a_num(vel_data_opt)
  B_opt = get_B_num(dtrq_data_opt)
  c_opt = get_c_num(acc_data_opt)
  D_opt = get_D_num(dtrq_data_opt, ddtrq_data_opt)
  vel_m_opt_2d = a_num + B_num @ alpha_s.value
  acc_m_opt_2d = c_num + D_num @ alpha_s.value
  vel_m_opt_2d = np.vstack((vel_m_opt_2d, vel_m_opt_2d[0][0]))
  acc_m_opt_2d = np.vstack((acc_m_opt_2d, acc_m_opt_2d[0][0]))
  vel_m_opt = np.zeros(data_len)
  acc_m_opt = np.zeros(data_len)
  for i in range(0, data_len):
    vel_m_opt[i]=vel_m_opt_2d[i][0]
    acc_m_opt[i]=acc_m_opt_2d[i][0]
  tau_m_opt = I_m*acc_m_opt+b_m*vel_m_opt-(1/r)*trq_data
  Em_true = integrate.simpson((1/(K_m*K_m))*tau_m_opt**2 + tau_m_opt*vel_m_opt, t)
  Em_per = (Em_true-Em_0)*100/Em_0
  print(f"Em_obj: {prob.value:.6f}; Em_true: {Em_true:.6f} ({Em_per:.1f}%)")

  Em_new = np.zeros(mc_no)
  vel_data_opt_rms = np.sqrt(np.mean(vel_data_opt**2))
  acc_data_opt_rms = np.sqrt(np.mean(acc_data_opt**2))
  trq_data_rms = np.sqrt(np.mean(trq_data**2))
  dtrq_data_opt_rms = np.sqrt(np.mean(dtrq_data_opt**2))
  ddtrq_data_opt_rms = np.sqrt(np.mean(ddtrq_data_opt**2))
  for j in range(0, mc_no):
    vel_data_new = vel_data_opt + vel_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    acc_data_new = acc_data_opt + acc_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    trq_data_new = trq_data + trq_data_rms*2*error*(np.random.rand(data_len)-0.5)
    dtrq_data_new = dtrq_data_opt + dtrq_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    ddtrq_data_new = ddtrq_data_opt + ddtrq_data_opt_rms*2*error*(np.random.rand(data_len-1)-0.5)
    a_new_1d = get_a_num(vel_data_new)
    B_new = get_B_num(dtrq_data_new)
    c_new_1d = get_c_num(acc_data_new)
    D_new = get_D_num(dtrq_data_new, ddtrq_data_new)
    a_new = a_new_1d[:, np.newaxis]
    c_new = c_new_1d[:, np.newaxis]

    vel_m_new_2d = a_new + B_new @ alpha_s.value
    acc_m_new_2d = c_new + D_new @ alpha_s.value
    vel_m_new_2d = np.vstack((vel_m_new_2d, vel_m_new_2d[0][0]))
    acc_m_new_2d = np.vstack((acc_m_new_2d, acc_m_new_2d[0][0]))
    vel_m_new = np.squeeze(vel_m_new_2d)
    acc_m_new = np.squeeze(acc_m_new_2d)
    
    tau_m_new = I_m*acc_m_new+b_m*vel_m_new-(1/r)*trq_data
    dq_m_new = vel_m_new
    Em_new[j] = integrate.simpson((1/(K_m*K_m))*tau_m_new**2 + tau_m_new*dq_m_new, t)

  Em_trial = np.mean(Em_new)
  Em_std = np.std(Em_new)

  Em_act_per = (Em_trial-Em_0)*100/Em_0
  Em_std_act = (Em_trial+2*Em_std-Em_0)*100/Em_0 - Em_act_per
  print(f"Em_trial: {Em_act_per:.1f}%; stdev: {Em_std_act:.1f}%")
  tidx = np.argmin(np.abs(TRQ_s))
  delta_s = delta_s - delta_s[tidx]
  axs.plot(delta_s, TRQ_s, 'o', label = 'nominal nonlinear')
  axs.legend()
  axs.set_xlabel('delta_s [rad]')
  axs.set_ylabel('tau_s [N-m]')
  plt.show()
  

  log_file_name = os.path.join('logs', 'sealog_' + cadence + '_' + joint + '_' + datetime.now().strftime("%m%d_%H%M%S") + '.txt')
  with open(log_file_name, 'a') as f:
    for idx in range(0, data_len-1):
      f.write(f'{alpha_s.value[idx][0]},\n')


