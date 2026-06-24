import numpy as np
import matplotlib.pyplot as plt
from melt_model import Model

# ---- user settings ----
masses = np.arange(1, 11, 1)   # total mass in Mars masses
angles = [0, 30, 45, 60, 90]
gamma = 0.1
vel = 1.5                     # vimp = vesc
entropy0 = 3160
Tmelt_model = "old_model"

def get_fractional_melt_volume(model, data):
    """
    User-defined output.

    Option 1: use model.du_melt after run_model().
    This is a simple grid-average placeholder.
    Replace this with a proper volume-weighted integral if needed.
    """
    return np.nanmean(model.du_melt)

    # Option 2, if you want a value from data:
    # return data["mantle melt mass fraction"]

# ---- run models ----
results = {}

for angle in angles:
    yvals = []

    for Mtotal in masses:
        print(f"Running Mtotal={Mtotal}, angle={angle}")

        m = Model(
            Mtotal=Mtotal,
            gamma=gamma,
            vel=vel,
            entropy0=entropy0,
            impact_angle=angle,
            Tmelt_model=Tmelt_model,
            use_tex=False,
        )

        data = m.run_model()
        yvals.append(get_fractional_melt_volume(m, data))

    results[angle] = np.array(yvals)

# ---- single combined plot ----
plt.figure(figsize=(6, 4))

for angle in angles:
    plt.plot(
        masses,
        results[angle],
        marker="o",
        label=f"{angle}°"
    )

plt.xlabel("Total mass (Mars mass)")
plt.ylabel("Fractional melt volume")
plt.title(f"gamma={gamma}, vimp=vesc")
plt.ylim(0, 1)
plt.legend(title="Impact angle")
plt.tight_layout()
plt.savefig("fig_fractional_melt_volume_all_angles_deng.png", dpi=300)
plt.show()
