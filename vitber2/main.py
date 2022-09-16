import ode
import animate
import numpy as np

y0 = [0, 0, 0, 0, 0.2, 0, 0, 3] #=[vx, vy, x, y, omega=0.2, theta, vs, s=3]
xs, ys = ode.solve_boat_dynamics(0, 40 , y0, mass_load= 0.001 * ode.m, fence=True)
red_xs = np.array([xs[i] for i in range(0, len(xs), 8)])
red_ys = np.array([ys[i] for i in range(0, len(ys), 8)])
print(red_ys)

animate.animate_deck_movement(red_xs, red_ys[:, 5], red_ys[:, 2], red_ys[:, 3], red_ys[:, 7], True)