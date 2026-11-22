"""
Generate textures/venus_clouds.jpg - what Venus looks like from space.

Venus is covered by thick sulfuric-acid clouds. In enhanced spacecraft images
(e.g. Mariner 10) it looks white to pale blue-grey, with long cloud streaks
that sweep from the equator toward the poles in a sideways "V", a warm tan
haze near the equator, and darker blue-grey toward the poles.
(venus.jpg is a radar map of the hidden surface, which can't be seen from outside.)

Run from the solar-system-sim folder:  python textures/make_venus_clouds.py
Needs numpy and Pillow (both come with Anaconda).
"""
import numpy as np
from PIL import Image

W, H = 2048, 1024
rng = np.random.default_rng(11)


def value_noise(u, v, fx, fy):
    """Smooth random noise at (u, v) in [0,1]; wraps around in u (longitude)"""
    grid = rng.random((fy + 2, fx))
    x = (u % 1.0) * fx
    y = np.clip(v, 0, 1) * fy
    x0 = np.floor(x).astype(int); y0 = np.floor(y).astype(int)
    tx = x - x0; ty = y - y0
    tx = tx * tx * (3 - 2 * tx); ty = ty * ty * (3 - 2 * ty)   # smoothstep
    x1 = (x0 + 1) % fx; x0 = x0 % fx; y1 = y0 + 1
    top = grid[y0, x0] * (1 - tx) + grid[y0, x1] * tx
    bot = grid[y1, x0] * (1 - tx) + grid[y1, x1] * tx
    return top * (1 - ty) + bot * ty


def fbm(u, v, fx, fy, octaves=4):
    """Several octaves of noise added together, normalised to roughly 0..1"""
    total, amp, norm = 0.0, 1.0, 0.0
    for _ in range(octaves):
        total = total + amp * value_noise(u, v, fx, fy)
        norm += amp
        amp *= 0.5; fx *= 2; fy *= 2
    return total / norm


u = (np.arange(W) + 0.5)[None, :] / W          # longitude 0..1
v = (np.arange(H) + 0.5)[:, None] / H          # 0 = north pole, 1 = south pole
lat = (0.5 - v) * 2                            # +1 north, -1 south

# Cloud streaks: noise stretched along longitude, bent poleward in a sideways V
bend = np.sign(lat) * np.abs(lat) ** 0.8 * (0.5 + 0.5 * np.cos(2 * np.pi * (u - 0.35)))
warp = 0.11 * bend + 0.05 * (fbm(u, v, 5, 4, 4) - 0.5) + 0.015 * (fbm(u, v, 16, 12, 3) - 0.5)
streaks = fbm(u, v + warp, 3, 48, 4)
streaks = np.clip((streaks - 0.5) * 2.2 + 0.5, 0, 1)

# Warm tan haze in patches around the equator
haze = fbm(u, v + 0.5 * warp, 5, 4, 5)
haze = np.clip((haze - 0.42) * 2.2, 0, 1) ** 1.3 * np.exp(-(lat / 0.45) ** 2)
haze *= 0.6 + 0.4 * streaks                    # haze shows through the lighter lanes

white = np.array([240, 240, 252], dtype=np.float32)   # (a bit blue: the sim's sunlight is slightly yellow)
streak_dark = np.array([192, 200, 222], dtype=np.float32)   # blue-grey cloud lanes
tan = np.array([212, 172, 128], dtype=np.float32)
polar = np.array([165, 176, 204], dtype=np.float32)

rgb = white + (streak_dark - white) * (1 - streaks)[..., None] * 0.8
rgb = rgb + (tan - rgb) * (haze * 0.8)[..., None]
rgb = rgb + (polar - rgb) * (np.abs(lat) ** 3 * 0.6)[..., None]
rgb = np.clip(rgb, 0, 255).astype(np.uint8)

Image.fromarray(rgb).save('textures/venus_clouds.jpg', quality=92)
print('saved textures/venus_clouds.jpg', rgb.reshape(-1, 3).mean(0).round())
