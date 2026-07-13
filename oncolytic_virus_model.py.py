import numpy as np
import matplotlib.pyplot as plt

HOURS = 72
DT = 0.1
steps = int(HOURS / DT)
time = np.linspace(0, HOURS, steps)

r = 0.04
Cmax = 1.0e7
k = 2.0e-7
p = 5.0e-7
d = 0.02

C0 = 1.0e6
V0 = 1.0e4

C = np.zeros(steps)
V = np.zeros(steps)
C[0] = C0
V[0] = V0

for i in range(steps - 1):
    dC = (r * C[i] * (1 - C[i] / Cmax) - k * V[i] * C[i]) * DT
    dV = (p * V[i] * C[i] - d * V[i]) * DT
    C[i + 1] = max(C[i] + dC, 0)
    V[i + 1] = max(V[i] + dV, 0)

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.set_xlabel("Time (hours)", fontsize=12)
ax1.set_ylabel("Cancer cells", color="#d62728", fontsize=12)
ax1.plot(time, C, color="#d62728", linewidth=2.2, label="Cancer cells")
ax1.tick_params(axis="y", labelcolor="#d62728")
ax1.set_xlim(0, HOURS)
ax1.set_ylim(bottom=0)

ax2 = ax1.twinx()
ax2.set_ylabel("Viral particles", color="#1f77b4", fontsize=12)
ax2.plot(time, V, color="#1f77b4", linewidth=2.2, label="Viral particles")
ax2.tick_params(axis="y", labelcolor="#1f77b4")
ax2.set_ylim(bottom=0)

plt.title("Oncolytic Virus Destroying Cancer Cells Over 72 Hours", fontsize=14, pad=15)
fig.tight_layout()
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=11)
plt.savefig("oncolytic_virus_graph.png", dpi=150)
plt.show()

print("Done! Graph saved as oncolytic_virus_graph.png")
print(f"Starting cancer cells: {C[0]:,.0f}")
print(f"Final cancer cells:    {C[-1]:,.0f}")
print(f"Peak viral particles:  {V.max():,.0f}")
