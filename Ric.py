import numpy as np
PI = np.pi
def Ric(x): 
  arg=0.5*x
  arg *= arg;
  return (1.-2.*arg)*np.exp(-arg)

def tri(x):
  if np.abs(x)>1: return 0
  return 1-np.abs(x)

def L2(x):
  if np.abs(x)>1.5: return 0
  if np.abs(x)>0.5: return 0.125*(3-2*np.abs(x))**2
  return 0.75-x**2

F0 = 0.02
Nw = 20
N = 1200 
dt = 1*0.2/Nw/F0
do = 2*np.pi/N/dt
o0 = 2*np.pi*F0
T0 = 2*np.pi/o0
Rm = dt*Nw
t1start = 3*T0 
t2start = 2*T0
M = N/Nw
print 'M = ',M
print 'dt,do=',dt,do,'\nomega_0, T_0 =', o0, T0, '\npulse 1,2 start time', t1start, t2start

fx = open("x.dat",'w') 
fX = open("X.dat",'w')
fm = open("M.dat",'w')

x1 = [Ric(o0*(dt*i-t1start)) for i in range(N)]
x2 = [Ric(o0*(dt*i-t2start)) for i in range(N)]

X1 = np.fft.fft(x1)
X2 = np.fft.fft(x2)
CXX = np.multiply(X1,X2)
cxx = np.fft.ifft(CXX)
for k in range(N):
  ok = k*do
  print >>fX, ok, np.abs(X1[k]), np.abs(X2[k]), np.abs(CXX[k])

A1m = np.zeros(M) 
P1m = np.zeros(M) 
P2m = np.zeros(M) 
A2m = np.zeros(M) 
for m in range(M):  #STFT
  x1m = [x1[n] for n in range(m*Nw,(m+1)*Nw)]
  x2m = [x2[n] for n in range(m*Nw,(m+1)*Nw)]
  cs1,cs2 = 0,0
  sn1,sn2 = 0,0
  for n in range(Nw):
    arg = dt*(n+m*Nw)*o0
    cs1+= x1m[n]*np.cos(arg)
    sn1+= x1m[n]*np.sin(arg)
    cs2+= x2m[n]*np.cos(arg)
    sn2+= x2m[n]*np.sin(arg)
  A1m[m] = np.sqrt(sn1*sn1+cs1*cs1)
  A2m[m] = np.sqrt(sn2*sn2+cs2*cs2)
  P1m[m] = np.arctan(sn1/cs1) if (A1m[m]>1e-4) else 0
  P2m[m] = np.arctan(sn2/cs2) if (A2m[m]>1e-4) else 0

#addP = 0
#for m in range(1,M-1):
  #P1m[m] +=addP*np.pi
  #if P1m[m+1]<P1m[m]: 
    #addP+=1
for m in range(M):  #plot STFT result in M.dat
  print >>fm, A1m[m], P1m[m]

ax1 = np.zeros(N)
for n in range(N): #Recover x1 to ax1. Compare recovered with original: plot  "x.dat" u 1:2, "" u 1:5
   mm0p = [n/Nw-1*(n/Nw>0),n/Nw,n/Nw+1*(n/Nw<M-1)]
   An,Pn = 0,0
   Am0p = [A1m[i] for i in mm0p]
   Pm0p = [P1m[i] for i in mm0p]
   for im in range(3):
       An += Am0p[im]*tri(1.0*n/Nw-mm0p[im]-0.5)
#       Pn += Pm0p[im]*tri(1.0*n/Nw-mm0p[im]-0.5) 
   ax1[n] = 2*An*np.cos(o0*n*dt-Pn)/Nw
   #ax1[n] = Pn

caxx = np.zeros(N)
for n in range(N): #Recover convolution fo x1 and x2 to to caxx. Compare recovered with original: plot  "x.dat" u 1:4, "" u 1:6
   An,Pn = 0,0
   for m1 in range(M):
    for m2 in (n/Nw-m1-1,n/Nw-m1,n/Nw-m1+1):
       if m2>0 and m2<M-2L:
        An += A1m[m1]*A2m[m2]*tri(1.0*n/Nw-m1-m2-1)
   caxx[n] = An/16.*np.cos(o0*n*dt+Pn)


for i in range(N):
   print>>fx, i*dt, x1[i], x2[i], cxx[i].real, ax1[i], caxx[i]

