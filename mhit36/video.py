#!/usr/bin/env python3
# 3D movie of the fiber in the periodic box, from output/fib_XXXXXXXX.dat (control points).
# Positions are wrapped into [0, lx)^3; where the fiber crosses a face of the box it is
# split and continues from the opposite face. Optionally shows the trajectory of the centre.
#
# usage: python3 video.py [--fps N] [--out fiber.mp4] [--dir output] [--lx 6.2831853]
#                           [--trail N] [--rotate DEG] [--elev 20] [--azim -60]
#   --trail   number of previous frames of the centre trajectory to draw (0 = none)
#   --rotate  rotation of the camera around the vertical axis over the whole movie (degrees)
#   output is .mp4 if ffmpeg is available, otherwise .gif
# needs numpy and matplotlib (e.g. ~/venvs/fibre/bin/python video.py)
import argparse, glob, math, os, re, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation

ap = argparse.ArgumentParser()
ap.add_argument('--fps', type=int, default=25)
ap.add_argument('--out', default='fiber.mp4')
ap.add_argument('--dir', default='output')
ap.add_argument('--lx', type=float, default=2.0*math.pi)
ap.add_argument('--trail', type=int, default=200)
ap.add_argument('--rotate', type=float, default=0.0)
ap.add_argument('--elev', type=float, default=20.0)
ap.add_argument('--azim', type=float, default=-60.0)
# parse_known_args: ignore the extra arguments added by Jupyter (VS Code interactive window)
args, _ = ap.parse_known_args()

files = sorted(glob.glob(os.path.join(args.dir, 'fib_[0-9]*.dat')))
if len(files) < 2:
    raise SystemExit('need at least two fib_XXXXXXXX.dat files in %s' % args.dir)
steps = [int(re.search(r'fib_(\d+)\.dat', f).group(1)) for f in files]
X = np.array([np.loadtxt(f, usecols=(0, 1, 2)) for f in files])   # (frame, node, xyz)
lx = args.lx
mid = (X.shape[1]+1)//2 - 1

# time step from input.inp (6th line), if available
dt = None
if os.path.exists('input.inp'):
    try:
        dt = float(open('input.inp').read().split('\n')[5].split()[0])
    except (IndexError, ValueError):
        dt = None

def wrapped_pieces(P):
    """Split the polyline P (nodes, 3) at the faces of the box, return pieces in [0,lx)^3."""
    # periodic image of each node; cut the fiber where it changes image
    img = np.floor(P/lx)
    out, piece = [], [P[0] - img[0]*lx]
    for k in range(1, len(P)):
        if np.any(img[k] != img[k-1]):
            # point where the segment leaves the box of node k-1
            a, b = P[k-1], P[k]
            t_cut = 1.0
            for d in range(3):
                if img[k, d] != img[k-1, d]:
                    face = max(img[k, d], img[k-1, d])*lx
                    t_cut = min(t_cut, (face - a[d])/(b[d] - a[d]))
            c = a + t_cut*(b - a)
            piece.append(c - img[k-1]*lx)
            out.append(np.array(piece))
            piece = [c - img[k]*lx]
        piece.append(P[k] - img[k]*lx)
    out.append(np.array(piece))
    return out

# centre trajectory, wrapped, with breaks (nan) where it jumps across the box
C = np.mod(X[:, mid, :], lx)
jump = np.any(np.abs(np.diff(C, axis=0)) > 0.5*lx, axis=1)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')
# box edges
for i in (0, lx):
    for j in (0, lx):
        ax.plot([0, lx], [i, i], [j, j], color='0.6', lw=0.8)
        ax.plot([i, i], [0, lx], [j, j], color='0.6', lw=0.8)
        ax.plot([i, i], [j, j], [0, lx], color='0.6', lw=0.8)
ax.set_xlim(0, lx); ax.set_ylim(0, lx); ax.set_zlim(0, lx)
ax.set_box_aspect((1, 1, 1))
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
ax.view_init(elev=args.elev, azim=args.azim)
npieces = 4   # a fiber shorter than the box crosses at most three faces
fib = [ax.plot([], [], [], '-', color='C3', lw=2.5)[0] for _ in range(npieces)]
trail, = ax.plot([], [], [], '-', color='C0', lw=0.8, alpha=0.7)
title = fig.suptitle('')
fig.tight_layout()

def draw(k):
    pieces = wrapped_pieces(X[k])
    for n, ln in enumerate(fib):
        if n < len(pieces):
            p = pieces[n]
            ln.set_data(p[:, 0], p[:, 1]); ln.set_3d_properties(p[:, 2])
        else:
            ln.set_data([], []); ln.set_3d_properties([])
    if args.trail > 0:
        k0 = max(0, k - args.trail)
        T = C[k0:k+1].copy()
        cut = np.where(jump[k0:k])[0]
        T = np.insert(T, cut + 1, np.nan, axis=0)
        trail.set_data(T[:, 0], T[:, 1]); trail.set_3d_properties(T[:, 2])
    if args.rotate:
        ax.view_init(elev=args.elev, azim=args.azim + args.rotate*k/max(1, len(files)-1))
    t = '   t = %.3f' % (steps[k]*dt) if dt else ''
    title.set_text('step %d%s' % (steps[k], t))
    return fib + [trail, title]

anim = animation.FuncAnimation(fig, draw, frames=len(files), interval=1000/args.fps)
out = args.out
if out.endswith('.mp4') and shutil.which('ffmpeg') is None:
    out = os.path.splitext(out)[0] + '.gif'
    print('ffmpeg not found, writing a gif instead')
if out.endswith('.mp4'):
    writer = animation.FFMpegWriter(fps=args.fps)
else:
    writer = animation.PillowWriter(fps=args.fps)
anim.save(out, writer=writer, dpi=110)
print('written %s: %d frames' % (out, len(files)))
