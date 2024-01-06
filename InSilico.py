'''
2023/12/12
created by Keishi Nishio (BGA23068)
'''

import math
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

glagh_limit=300

class GFPExpressionCalculator:
    def __init__(self, alpha_L, alpha_C, alpha_L_star, alpha_G, gamma_L, gamma_C, gamma_G, beta_A, beta_C, beta_L, n_A, n_C, n_L):
        self.alpha_L = alpha_L
        self.alpha_C = alpha_C
        self.alpha_L_star = alpha_L_star
        self.alpha_G = alpha_G
        self.gamma_L = gamma_L
        self.gamma_C = gamma_C
        self.gamma_G = gamma_G
        self.beta_A = beta_A
        self.beta_C = beta_C
        self.beta_L = beta_L
        self.n_A = n_A
        self.n_C = n_C
        self.n_L = n_L

    def calculate_GFP(self, A):
        def steady_state(variables):
            L_hat, C_hat, L_star_hat, G_hat = variables
            dLdt = self.alpha_L * A**self.n_A / (A**self.n_A + self.beta_A**self.n_A) - self.gamma_L * L_hat
            dCdt = self.alpha_C * A**self.n_A / (A**self.n_A + self.beta_A**self.n_A) - self.gamma_C * C_hat
            dL_stardt = self.alpha_L_star * self.beta_C**self.n_C / (self.beta_C**self.n_C + C_hat**self.n_C) - self.gamma_L * L_star_hat
            dGFPdt = self.alpha_G * self.beta_L**self.n_L / (self.beta_L**self.n_L + (L_hat + L_star_hat)**self.n_L) - self.gamma_G * G_hat
            return [dLdt, dCdt, dL_stardt, dGFPdt]

        initial_guesses = [0.1, 0.1, 0.1, 0.1]
        return fsolve(steady_state, initial_guesses)

class AHLConcentrationCalculator:
    def __init__(self, lambda_val):
        self.lambda_val = lambda_val
        self.points = []

    def add_point(self, x, y, A_0):
        self.points.append((x, y, A_0))

    def distance_to(self, x, y, other_x, other_y):
        return math.sqrt((x - other_x) ** 2 + (y - other_y) ** 2)

    def calc_sum_AHL_concentration(self, x, y):
        total_concentration = 0
        for point in self.points:
            distance = self.distance_to(x, y, point[0], point[1])
            total_concentration += point[2] * math.exp(-distance / self.lambda_val)
        return total_concentration

    def calc_concentration_grid(self):
        concentration_grid = {}
        for x in np.arange(-1*glagh_limit, glagh_limit+1, 5):
            for y in np.arange(-1*glagh_limit, glagh_limit+1, 5):
                concentration = self.calc_sum_AHL_concentration(x, y)
                concentration_grid[(round(x, 1), round(y, 1))] = concentration
        return concentration_grid

class GFPHighConcentrationLocator:
    def __init__(self, gfp_calculator, ahl_calculator):
        self.gfp_calculator = gfp_calculator
        self.ahl_calculator = ahl_calculator

    def find_high_concentration_positions(self, threshold):
        GFP_high_positions = []
        concentration_grid = self.ahl_calculator.calc_concentration_grid()
        for position, A in concentration_grid.items():
            GFP = self.gfp_calculator.calculate_GFP(A)
            if GFP[3] >= threshold:
                GFP_high_positions.append(position)
        return GFP_high_positions

# Parameters
alpha_L= 1
alpha_C= 1
alpha_L_star= 1
alpha_G= 1
gamma_L= 0.1
gamma_C= 0.1
gamma_G= 0.1
beta_A= 0.01
beta_C= 0.01
beta_L= 1
n_A= 2
n_C= 1
n_L=2
lambda_val = 10

# Create instances of the classes
gfp_calculator = GFPExpressionCalculator(alpha_L, alpha_C, alpha_L_star, alpha_G, gamma_L, gamma_C, gamma_G, beta_A, beta_C, beta_L, n_A, n_C, n_L)
ahl_calculator = AHLConcentrationCalculator(lambda_val)


########################################################
# Add points with specific initial concentrations


# pattern A
# ahl_calculator.add_point(0, 0, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
# ahl_calculator.add_point(60, 0, 5)   # 位置 (10, 0) に A_0 = 1 のポイントを追加

# pattern B
# ahl_calculator.add_point(0, 0, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
# ahl_calculator.add_point(60, 0, 20)   # 位置 (10, 0) に A_0 = 1 のポイントを追加

# pattern C
# ahl_calculator.add_point(0, 0, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
# ahl_calculator.add_point(60, 0, 5)   # 位置 (10, 0) に A_0 = 1 のポイントを追加
# ahl_calculator.add_point(30, -30, 5)   # 位置 (10, 0) に A_0 = 1 のポイントを追加

# pattern D
ahl_calculator.add_point(-60, -60, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
ahl_calculator.add_point(-60, 60, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
ahl_calculator.add_point(60, -60, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
ahl_calculator.add_point(60, 60, 5)   # 位置 (0, 0) に A_0 = 10 のポイントを追加
#########################################################



# Create a locator instance and find positions with high GFP concentration
gfp_locator = GFPHighConcentrationLocator(gfp_calculator, ahl_calculator)
GFP_high_positions = gfp_locator.find_high_concentration_positions(2)

# Plotting
# ... [Previous code]

# Plotting
plt.figure(figsize=(10, 10))

# Plotting positions with high GFP concentration
for pos in GFP_high_positions:
    plt.plot(pos[0], pos[1], 'go')  # Green circles for positions with high GFP concentration

# Plotting the added points
for point in ahl_calculator.points:
    x, y, A_0 = point
    print(x)
    print(y)
    print(A_0)
    plt.plot(x, y, 'ro', markersize=A_0*2)  # Red dot for added points, size proportional to A_0

plt.xlim(-1*glagh_limit, glagh_limit)
plt.ylim(-1*glagh_limit, glagh_limit)
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Positions with GFP Concentration >= 2')
plt.grid(True)

# plt.savefig("patternA")
# plt.savefig("patternB")
# plt.savefig("patternC")
plt.savefig("patternD")
plt.show()

