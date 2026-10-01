#!/usr/bin/env python3
# Natural-frequency test of the fiber (fib_coupling=.false., pert_amp>0, fib_log=1).
# Measures the first bending frequency from the zero crossings of the displacement
# of the fiber centre and compares it with the free-free Euler-Bernoulli value.
#
# usage: python3 fiber_freq.py [dt] [component] [log file]
#   dt        time step of the run (default 2.5e-4)
#   component 0, 1 or 2 = x, y or z displacement of the centre (default 1, as pert_dir=(0,1,0))
import sys, math

dt   = float(sys.argv[1]) if len(sys.argv) > 1 else 2.5e-4
comp = int(sys.argv[2])   if len(sys.argv) > 2 else 1
fname = sys.argv[3]       if len(sys.argv) > 3 else 'output/fiber_log.dat'

# fiber properties (keep in sync with fib_param.f90)
L, d, rhof = 1.0, 0.025, 10.0
E, G = 218.0e3, 83.846e3

r = 0.5*d
A = math.pi*r**2
I = math.pi/4.0*r**4
bL = 4.730040744862704
omega = lambda EI: bL**2*math.sqrt(EI/(rhof*A*L**4))

# log columns: step, NR iterations, length, centre x y z, centre u v w
rows = [list(map(float, l.split())) for l in open(fname) if l.strip()]
t = [row[0]*dt for row in rows]
x = [row[3+comp] for row in rows]
xm = sum(x)/len(x)
x = [v - xm for v in x]
tz = [t[i] - x[i]*(t[i+1]-t[i])/(x[i+1]-x[i]) for i in range(len(x)-1) if x[i]*x[i+1] < 0]
if len(tz) < 3:
    sys.exit('not enough oscillations in %s' % fname)
T = 2.0*(tz[-1]-tz[0])/(len(tz)-1)
w = 2.0*math.pi/T

print('numerical  omega1 = %.4f   (%d zero crossings, %d steps)' % (w, len(tz), len(rows)))
print('analytical omega1 = %.4f   (E*I)   error %+.2f%%' % (omega(E*I), 100*(w/omega(E*I)-1)))
print('analytical omega1 = %.4f   (G*J)   error %+.2f%%' % (omega(G*2*I), 100*(w/omega(G*2*I)-1)))
print('max length drift  = %.2e,  max Newton iterations = %d'
      % (max(abs(row[2]-L) for row in rows), int(max(row[1] for row in rows))))
