import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)
pd.set_option('display.expand_frame_repr', False)

B, D = 0.1, 0.01
tau_b, tau_d = 0.3, 0.1
m = n = 5
xi = eta = 0.01

default_tau_a = 1.0

def f_act(s_im1, s_i, s_ip1):
    Bi = B * (1.0 + s_ip1) / (1.0 + xi * s_ip1) * (1.0 + eta * s_i)
    return 1.0 / (1.0 + (Bi / max(s_im1, 1e-16))**m)

def g_sup(s_im1, s_i, s_ip1):
    Di = D * (1.0 + s_im1) / (1.0 + xi * s_im1) * (1.0 + eta * s_i)
    return 1.0 / (1.0 + (Di / max(s_ip1, 1e-16))**n)

def rhs(t, y, tau_a):
    s1, s2, s3, s4, a = y
    dead = 1e-4
    if s1 < dead: s1 = 0.0
    if s2 < dead: s2 = 0.0
    if s3 < dead: s3 = 0.0
    if s4 < dead: s4 = 0.0

    s0 = a
    s5 = 0.0

    F1 = (1.0/tau_b)*f_act(s0, s1, s2)*s1 - (1.0/tau_d)*g_sup(s0, s1, s2)*s1
    F2 = (1.0/tau_b)*f_act(s1, s2, s3)*s2 - (1.0/tau_d)*g_sup(s1, s2, s3)*s2
    F3 = (1.0/tau_b)*f_act(s2, s3, s4)*s3 - (1.0/tau_d)*g_sup(s2, s3, s4)*s3
    F4 = (1.0/tau_b)*f_act(s3, s4, s5)*s4 - (1.0/tau_d)*g_sup(s3, s4, s5)*s4

    elim = (1.0/tau_a) * (1.0 / (1.0 + (1.0/max(s1,1e-16))**2)) * a
    da = -elim

    return np.array([F1, F2, F3, F4, da], float)

def rk4_step(t, y, h, tau_a):
    k1 = rhs(t, y, tau_a)
    k2 = rhs(t + 0.5*h, y + 0.5*h*k1, tau_a)
    k3 = rhs(t + 0.5*h, y + 0.5*h*k2, tau_a)
    k4 = rhs(t + h, y + h*k3, tau_a)
    yn = y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    yn[:4] = np.maximum(yn[:4], 0.0)
    mask = yn[:4] < 1e-4
    yn[:4][mask] = 0.0
    yn[4] = max(yn[4], 0.0)
    return yn

def simulate(a0, T=12.0, dt=1e-3, tau_a=default_tau_a):
    y = np.array([1e-3, 1e-3, 1e-3, 1e-3, a0], float)
    steps = int(np.ceil(T/dt))
    out = np.empty((steps+1, 6), float)
    t = 0.0
    out[0] = [t, *y]
    for k in range(1, steps+1):
        y = rk4_step(t, y, dt, tau_a)
        t = k*dt
        out[k] = [t, *y]
    return out

def my_plot(traj, title, fname, ylim=(-4.5, 4.5)):
    t, s1, s2, s3, s4, a = traj.T
    eps = 1e-12
    L1 = np.log10(np.maximum(s1, eps))
    L2 = np.log10(np.maximum(s2, eps))
    L3 = np.log10(np.maximum(s3, eps))
    L4 = np.log10(np.maximum(s4, eps))
    La = np.log10(np.maximum(a,  eps))

    plt.figure(figsize=(7,5))
    plt.plot(t, L1, label="Ab-1")
    plt.plot(t, L2, label="Ab-2")
    plt.plot(t, L3, label="Ab-3")
    plt.plot(t, L4, label="Ab-4")
    plt.plot(t, La, linestyle='--', label="a (antigen)")
    plt.xlabel("Time")
    plt.ylabel("log₁₀(concentration)")
    plt.xlim(0, 12)
    plt.ylim(*ylim)
    plt.yticks([-4,-3,-2,-1,0,1,2])
    plt.grid(True)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close()

def summarize(traj):
    t, s1, s2, s3, s4, a = traj.T
    def peak_and_time(x):
        i = int(np.argmax(x))
        return float(np.max(x)), float(t[i])
    peaks = {}
    for name, s in zip(["Ab1","Ab2","Ab3","Ab4"], [s1,s2,s3,s4]):
        peaks[f"peak_{name}"], peaks[f"t_peak_{name}"] = peak_and_time(s)
    peaks["a_final"] = float(a[-1])
    return peaks

if __name__ == "__main__":
    runs = [
        {"a0": 1.0,   "tau_a": 10.0, "name": "a0_1e0"},   # LZT
        {"a0": 100.0, "tau_a": 1.0,  "name": "a0_1e2"},   # NRP
        {"a0": 1000.0,"tau_a": 0.5,  "name": "a0_1e3"},   # HZT
    ]

    summary_data = {}

    for r in runs:
        a0, tauA, name = r["a0"], r["tau_a"], r["name"]
        traj = simulate(a0=a0, T=12.0, dt=1e-3, tau_a=tauA)
        summary_data[name] = summarize(traj)
        my_plot(traj, f"a(0)={a0:g}, τ_a={tauA:g}", f"{name}_plot.png")
        print(f"[ok] saved: {name}_plot.png")

    df = pd.DataFrame(summary_data).T
    df = df[["peak_Ab1","t_peak_Ab1","peak_Ab2","t_peak_Ab2",
             "peak_Ab3","t_peak_Ab3","peak_Ab4","t_peak_Ab4","a_final"]]

    df.index = [
        r"$a(0)=10^{0}$",
        r"$a(0)=10^{2}$",
        r"$a(0)=10^{3}$"
    ]

    df.columns = [
        "Peak Ab1", "t_peak(Ab1)",
        "Peak Ab2", "t_peak(Ab2)",
        "Peak Ab3", "t_peak(Ab3)",
        "Peak Ab4", "t_peak(Ab4)",
        "a_final"
    ]

    print("\n===== Summary of immune responses =====")
    print(df.round(4).to_string())