#!/usr/bin/env python3
# Movie of the fiber from output/fib_XXXXXXXX.dat (control points, written every fib_out steps).
# Three orthogonal views of the fiber. Small displacements (e.g. the natural-frequency
# test) are magnified with respect to the first frame so that the bending is visible.
#
# usage: python3 fiber_movie.py [--scale S] [--fps N] [--out fiber.mp4] [--dir output]
#   --scale  magnification of the displacement from the first frame
#            (default: automatic, max displacement shown as 10% of the fiber length;
#             use --scale 1 to see the true shape)
#   output is .mp4 if ffmpeg is available, otherwise .gif
# needs numpy and matplotlib (e.g. ~/venvs/fibre/bin/python fiber_movie.py)
import argparse, glob, os, re, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation

ap = argparse.ArgumentParser()
ap.add_argument('--scale', type=float, default=None)
ap.add_argument('--fps', type=int, default=25)
ap.add_argument('--out', default='fiber.mp4')
ap.add_argument('--dir', default='output')
# parse_known_args: ignore the extra arguments added by Jupyter (VS Code interactive window)
args, _ = ap.parse_known_args()

files = sorted(glob.glob(os.path.join(args.dir, 'fib_[0-9]*.dat')))
if len(files) < 2:
    raise SystemExit('need at least two fib_XXXXXXXX.dat files in %s' % args.dir)
steps = [int(re.search(r'fib_(\d+)\.dat', f).group(1)) for f in files]
X = np.array([np.loadtxt(f, usecols=(0, 1, 2)) for f in files])   # (frame, node, xyz)

# time step from input.inp (6th line), if available
dt = None
if os.path.exists('input.inp'):
    try:
        dt = float(open('input.inp').read().split('\n')[5].split()[0])
    except (IndexError, ValueError):
        dt = None

# magnified displacement from the first frame
X0 = X[0]
L = np.sum(np.linalg.norm(np.diff(X0, axis=0), axis=1))
dmax = np.abs(X - X0).max()
scale = args.scale if args.scale is not None else (0.1*L/dmax if dmax > 0 else 1.0)
Y = X0 + scale*(X - X0)

views = [(0, 1, 'x', 'y'), (0, 2, 'x', 'z'), (1, 2, 'y', 'z')]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
lines = []
cen = 0.5*(Y.min(axis=(0, 1)) + Y.max(axis=(0, 1)))
half = 0.55*max((Y.max(axis=(0, 1)) - Y.min(axis=(0, 1))).max(), 1e-12)
for ax, (a, b, la, lb) in zip(axes, views):
    ax.plot(X0[:, a], X0[:, b], color='0.75', lw=1, ls='--', label='t = 0')
    ln, = ax.plot([], [], '-o', color='C0', lw=2, ms=3)
    lines.append(ln)
    ax.set_xlim(cen[a]-half, cen[a]+half)
    ax.set_ylim(cen[b]-half, cen[b]+half)
    ax.set_aspect('equal')
    ax.set_xlabel(la)
    ax.set_ylabel(lb)
title = fig.suptitle('')
fig.tight_layout(rect=(0, 0, 1, 0.93))

def draw(k):
    for ln, (a, b, _, _) in zip(lines, views):
        ln.set_data(Y[k, :, a], Y[k, :, b])
    t = '   t = %.4f' % (steps[k]*dt) if dt else ''
    title.set_text('step %d%s   (displacement x %.3g)' % (steps[k], t, scale))
    return lines + [title]

anim = animation.FuncAnimation(fig, draw, frames=len(files), interval=1000/args.fps, blit=False)
out = args.out
if out.endswith('.mp4') and shutil.which('ffmpeg') is None:
    out = os.path.splitext(out)[0] + '.gif'
    print('ffmpeg not found, writing a gif instead')
if out.endswith('.mp4'):
    writer = animation.FFMpegWriter(fps=args.fps)
else:
    writer = animation.PillowWriter(fps=args.fps)
anim.save(out, writer=writer, dpi=110)
print('written %s: %d frames, displacement magnified x %.3g' % (out, len(files), scale))
