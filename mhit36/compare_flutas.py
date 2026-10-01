#!/usr/bin/env python3
# Sedimenting fiber: MHIT36 vs FluTAS.
# FluTAS writes data/pos.txt every 100 steps (Update_Pos.f90):
#   time, |atan((z_end - z_c)/(y_end - y_c))|/pi, 2*sqrt(v_end^2 + w_end^2)
# with c the centre at the previous step and "end" the last control point.
# The same quantities are computed here from the MHIT36 files output/fib_XXXXXXXX.dat
# (the centre is taken at the same step, which makes no visible difference).
#
# usage: python3 compare_flutas.py [path/to/flutas/data/pos.txt] [--dt 5e-4] [--dir output]
import argparse, glob, math, os, re

ap = argparse.ArgumentParser()
ap.add_argument('pos', nargs='?', default='../flutas/src/data/pos.txt', help='FluTAS data/pos.txt')
ap.add_argument('--dt', type=float, default=5.0e-4)
ap.add_argument('--dir', default='output')
ap.add_argument('--png', default='compare_flutas.png')
# parse_known_args: ignore the extra arguments added by Jupyter (VS Code interactive window)
args, _ = ap.parse_known_args()
if args.pos and not os.path.exists(args.pos):
    print('FluTAS file %s not found: MHIT36 only' % args.pos)
    args.pos = None

# MHIT36
mh = []
for f in sorted(glob.glob(os.path.join(args.dir, 'fib_[0-9]*.dat'))):
    step = int(re.search(r'fib_(\d+)\.dat', f).group(1))
    rows = [list(map(float, l.split())) for l in open(f) if l.strip() and not l.startswith('#')]
    c, e = rows[(len(rows)+1)//2 - 1], rows[-1]          # centre and last control point
    ang = abs(math.atan((e[2]-c[2])/(e[1]-c[1]))/math.pi) if e[1] != c[1] else 0.5
    spd = 2.0*math.sqrt(e[8]**2 + e[9]**2)                # columns: x y z q1-q4 dxdt dydt dzdt
    mh.append((step*args.dt, ang, spd))
if not mh:
    raise SystemExit('no fib_XXXXXXXX.dat files in %s' % args.dir)

# FluTAS
fl = []
if args.pos:
    fl = [tuple(map(float, l.split()[:3])) for l in open(args.pos) if l.strip()]

print('%10s | %10s %10s | %10s %10s' % ('time', 'angle MH', 'speed MH', 'angle FL', 'speed FL'))
for t, a, s in mh[::max(1, len(mh)//25)]:
    if fl:
        k = min(range(len(fl)), key=lambda i: abs(fl[i][0]-t))
        print('%10.4f | %10.5f %10.5f | %10.5f %10.5f' % (t, a, s, fl[k][1], fl[k][2]))
    else:
        print('%10.4f | %10.5f %10.5f |' % (t, a, s))

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit('matplotlib not available: no plot')
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot([m[0] for m in mh], [m[2] for m in mh], label='MHIT36')
ax[1].plot([m[0] for m in mh], [m[1] for m in mh], label='MHIT36')
if fl:
    ax[0].plot([f[0] for f in fl], [f[2] for f in fl], '--', label='FluTAS')
    ax[1].plot([f[0] for f in fl], [f[1] for f in fl], '--', label='FluTAS')
ax[0].set_xlabel('t'); ax[0].set_ylabel('2 |(v,w)| at the fiber end')
ax[1].set_xlabel('t'); ax[1].set_ylabel('|angle in y-z plane| / pi')
for a in ax:
    a.legend(); a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(args.png, dpi=120)
print('plot written to', args.png)
