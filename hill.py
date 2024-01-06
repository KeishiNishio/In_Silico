'''
2023/12/12
created by Keishi Nishio (BGA23068)
'''

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class DecayModel:
    def __init__(self, A0, lambda_decay):
        self.A0 = A0
        self.lambda_decay = lambda_decay

    def calculate_A(self, d):
        return self.A0 * np.exp(-d / self.lambda_decay)


class SteadyStateModel:
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

    def calculate_steady_states(self, A):
        def steady_state(variables):
            L_hat, C_hat, L_star_hat, G_hat = variables
            dLdt = self.alpha_L * A**self.n_A / (A**self.n_A + self.beta_A**self.n_A) - self.gamma_L * L_hat
            dCdt = self.alpha_C * A**self.n_A / (A**self.n_A + self.beta_A**self.n_A) - self.gamma_C * C_hat
            dL_stardt = self.alpha_L_star * self.beta_C**self.n_C / (self.beta_C**self.n_C + C_hat**self.n_C) - self.gamma_L * L_star_hat
            dGdt = self.alpha_G * self.beta_L**self.n_L / (self.beta_L**self.n_L + (L_hat + L_star_hat)**self.n_L) - self.gamma_G * G_hat
            return [dLdt, dCdt, dL_stardt, dGdt]
        
        initial_guesses = [0.1, 0.1, 0.1, 0.1]
        return fsolve(steady_state, initial_guesses)


class Simulation:
    def __init__(self, decay_model, steady_state_model):
        self.decay_model = decay_model
        self.steady_state_model = steady_state_model

    def run_simulation(self, d_values):
        A_values = self.decay_model.calculate_A(d_values)
        steady_states = np.array([self.steady_state_model.calculate_steady_states(A) for A in A_values])
        return d_values, A_values, steady_states

    def plot_results(self, d_values, A_values, steady_states):
        L_hat_values, C_hat_values, L_star_hat_values, G_hat_values = steady_states.T

        plt.figure(figsize=(14, 8))

        plt.subplot(2, 2, 1)
        plt.plot(A_values, L_hat_values, label='L_hat')
        plt.xlim(0,0.025)
        plt.title('L_hat over A')
        plt.xlabel('A')
        plt.ylabel('L_hat')
        plt.grid(True)
        plt.legend()

        plt.subplot(2, 2, 2)
        plt.plot(A_values, C_hat_values, label='C_hat')
        plt.xlim(0,0.025)
        plt.title('C_hat over A')
        plt.xlabel('A')
        plt.ylabel('C_hat')
        plt.grid(True)
        plt.legend()

        plt.subplot(2, 2, 3)
        plt.plot(A_values, L_star_hat_values, label='L_star_hat')
        plt.xlim(0,0.025)
        plt.title('L_star_hat over A')
        plt.xlabel('A')
        plt.ylabel('L_star_hat')
        plt.grid(True)
        plt.legend()

        plt.subplot(2, 2, 4)
        plt.plot(A_values, G_hat_values, label='G_hat')
        plt.xlim(0,0.025)
        plt.title('G_hat over A')
        plt.xlabel('A')
        plt.ylabel('G_hat')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.savefig("hill.png")
        plt.show()
        

# 使用例
decay_model = DecayModel(A0=10, lambda_decay=10)
steady_state_model = SteadyStateModel(alpha_L=1, alpha_C=1, alpha_L_star=1, alpha_G=1, gamma_L=0.1, gamma_C=0.1, gamma_G=0.1, beta_A=0.01, beta_C=0.01, beta_L=1, n_A=2, n_C=1, n_L=2)
simulation = Simulation(decay_model, steady_state_model)

d_values, A_values, steady_states = simulation.run_simulation(np.linspace(0, 500, 10000))
simulation.plot_results(d_values, A_values, steady_states)
