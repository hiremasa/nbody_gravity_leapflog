import numpy as np

def generate_int_state():
  #初期条件の生成
  count=0
  n_xyz=[] #座標 r<1
  n_vel=[] #速度 |v|<1/8
  N=128 #粒子の数
  np.random.seed(0)
  while count<N:
    x=np.random.uniform(-1, 1)
    y=np.random.uniform(-1, 1)
    z=np.random.uniform(-1, 1)
    vx=np.random.uniform(-1/8, 1/8)
    vy=np.random.uniform(-1/8, 1/8)
    vz=np.random.uniform(-1/8, 1/8)
    if x**2 + y**2 + z**2<=1.0 and vx**2 + vy**2 + vz**2<=1/64:
      count+=1
      n_xyz.append([x,y,z])
      n_vel.append([vx,vy,vz])
  return n_xyz, n_vel

def calc_E_and_L(n_xyz, n_vel):
  sum_E=0
  sum_L=0
  eps=1e-2
  for i in range(N):
    sum_E+=0.5 * 1/N * sum(v**2 for v in n_vel[i]) +  0.5* 1/N * (1/N / eps)
    sum_L+=1/N * np.array([n_xyz[i][1]*n_vel[i][2] - n_xyz[i][2]*n_vel[i][1], n_xyz[i][2]*n_vel[i][0] - n_xyz[i][0]*n_vel[i][2],n_xyz[i][0]*n_vel[i][1] - n_xyz[i][1]*n_vel[i][0]])
    for j in range(N):
      dx=n_xyz[j][0]-n_xyz[i][0]
      dy=n_xyz[j][1]-n_xyz[i][1]
      dz=n_xyz[j][2]-n_xyz[i][2]
      r2=(dx**2 + dy**2 + dz**2+ eps**2)**(1/2)
      sum_E+= 0.5 * 1/N * (-1/N / r2) 
  print(sum_E)
  print(sum_L)


def calc_acc(n_xyz):
  n_acc=[]
  eps=1e-2
  for i in range(N):
    a_x = 0
    a_y = 0
    a_z = 0
    for j in range(N):
      dx=n_xyz[j][0]-n_xyz[i][0]
      dy=n_xyz[j][1]-n_xyz[i][1]
      dz=n_xyz[j][2]-n_xyz[i][2]
      r2=(dx**2 + dy**2 + dz**2+ eps**2)**(3/2)
      a_x+=(1/N) * dx / r2
      a_y+=(1/N) * dy / r2
      a_z+=(1/N) * dz / r2
    n_acc.append([a_x, a_y, a_z])
  return n_acc

N_xyz, N_vel = generate_int_state()
E0, _ = calc_E_and_L(N_xyz, N_vel)
delta_t=1/256
t=0
while t<1.0:
  new_xyz=[]
  new_vel=[]
  temp_v_12=[]
  N_acc=calc_acc(N_xyz)
  for i in range(N):
    x_new=N_xyz[i][0] + N_vel[i][0]*delta_t + 0.5 * N_acc[i][0]*delta_t**2
    y_new=N_xyz[i][1] + N_vel[i][1]*delta_t + 0.5 * N_acc[i][1]*delta_t**2
    z_new=N_xyz[i][2] + N_vel[i][2]*delta_t + 0.5 * N_acc[i][2]*delta_t**2

    vx_12=N_vel[i][0] + 0.5*N_acc[i][0]*delta_t
    vy_12=N_vel[i][1] + 0.5*N_acc[i][1]*delta_t
    vz_12=N_vel[i][2] + 0.5*N_acc[i][2]*delta_t

    new_xyz.append([x_new, y_new, z_new])
    temp_v_12.append([vx_12, vy_12, vz_12])

  new_acc=calc_acc(new_xyz)
  for i in range(N):
    vx_new=temp_v_12[i][0] + 0.5 * new_acc[i][0]*delta_t
    vy_new=temp_v_12[i][1] + 0.5 * new_acc[i][1]*delta_t
    vz_new=temp_v_12[i][2] + 0.5 * new_acc[i][2]*delta_t
    new_vel.append([vx_new, vy_new, vz_new])

  n_xyz=new_xyz
  n_vel=new_vel
  t+=delta_t

E, _ = calc_E_and_L(new_xyz, new_vel)

print(np.abs((E-E_0)/E_0))
