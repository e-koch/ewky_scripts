
import matplotlib.pyplot as plt
from itertools import cycle

'''
Some simple code for parsing the output abundances from Rate13 models.
There's also a simple plotting function for comparing the time evolution of
species abundances.
'''


def parse_output(filename):
    '''
    Parse the output files from Rate13.
    '''

    species_dict = {}

    with open(filename) as f:
        lines = f.readlines()

    # Now we'll parse through each set of abundances
    for i, line in enumerate(lines):
        splits = lines[i].split()
        if len(splits) == 0:
            continue

        if splits[0] == "TIME":
            keys = splits

            for key in keys:
                if key in species_dict.keys():
                    continue
                species_dict[key] = []

            continue

        for key, val in zip(keys, splits):
            if 'E' not in val:
                val = "E-".join(val.split("-"))
                if val[0] == "E":
                    val = val[2:]
            species_dict[key].append(float(val))

    species_dict["TIME"] = species_dict["TIME"][:80]

    return species_dict


def plot_species(species_dict, ax,
                 species=["C+", "C", "CO", "C2H", "HC3N",
                          "HC5N", "CH3OH", "e-"],
                 legend=True, loc='best', ylabel=True, xlabel=True,
                 title=None):
    '''
    Plot fraction of chemical species vs. time.
    '''

    linestyles = cycle(["-", "--", "-.", ":"])

    for spec, ls in zip(species, linestyles):
        ax.loglog(species_dict["TIME"], species_dict[spec], label=spec,
                  linestyle=ls)

    if legend:
        ax.legend(frameon=True, loc=loc)

    if xlabel:
        ax.set_xlabel("Age (yr)")
    if ylabel:
        ax.set_ylabel(r"Fraction of H$_2$")
    ax.grid()

    if title is not None:
        ax.set_title(title)


if __name__ == '__main__':

    import seaborn as sb
    sb.set_palette("colorblind", n_colors=10)

    # Make some plots out the output abundances in outputs

    default = parse_output("outputs/dc_0.out")
    freeze_out = parse_output("outputs/dc_1_freezeout.out")
    carbonrich = parse_output("outputs/dc_2_carbonrich.out")
    metalrich = parse_output("outputs/dc_3_metalrich.out")
    densecore = parse_output("outputs/dc_4_densecore.out")
    diffusecloud = parse_output("outputs/dc_5_diffusecloud.out")
    diffuseism = parse_output("outputs/dc_6_diffuseism.out")

    # Default comparisons
    plot_species(default, plt.subplot(111))
    plt.savefig("outputs/default_run.pdf")
    plt.close()

    # Effect of freeze-out
    fig, axes = plt.subplots(1, 2, sharey=True)
    plot_species(default, axes[0], title="Default")
    plot_species(freeze_out, axes[1], legend=False, ylabel=False,
                 title="Freeze-out")
    sb.despine()
    plt.savefig("outputs/freeze_out_affect.pdf")
    plt.close()

    # Vs. carbon-rich
    fig, axes = plt.subplots(1, 2, sharey=True)
    plot_species(default, axes[0], title="O-rich")
    plot_species(carbonrich, axes[1], legend=False, ylabel=False,
                 title="C-rich")
    sb.despine()
    plt.savefig("outputs/c_vs_o_rich.pdf")
    plt.close()

    # Vs. high-metallicity
    fig, axes = plt.subplots(1, 2, sharey=True)
    plot_species(default, axes[0], title="Default",
                 species=["CO", "SO", "NO", "SiO", "SO2"])
    plot_species(metalrich, axes[1], legend=False, ylabel=False,
                 title="Metal-rich", species=["CO", "SO", "NO", "SiO", "SO2"])
    sb.despine()
    plt.savefig("outputs/metal_rich.pdf")
    plt.close()

    # Comparing environments
    fig, axes = plt.subplots(2, 2, sharey=False)
    plot_species(densecore, axes[0, 0],
                 title=r"Dense Core (n=10$^5$ cm$^{-3}$, T=5 K, A$_{\rm v}$=15)",
                 species=["CO", "C", "C+", "e-"], xlabel=False)
    plot_species(default, axes[0, 1], legend=False, ylabel=False, xlabel=False,
                 title=r"Dense Cloud (n=10$^4$ cm$^{-3}$, T=10 K, A$_{\rm v}$=10)",
                 species=["CO", "C", "C+", "e-"])
    plot_species(diffusecloud, axes[1, 0], legend=False, ylabel=True,
                 title=r"Diffuse Cloud (n=10$^3$ cm$^{-3}$, T=50 K, A$_{\rm v}$=5)",
                 species=["CO", "C", "C+", "e-"])
    plot_species(diffuseism, axes[1, 1], legend=False, ylabel=False,
                 title=r"CNM (n=10$^2$ cm$^{-3}$, T=100 K, A$_{\rm v}$=2)",
                 species=["CO", "C", "C+", "e-"])
    axes[0, 0].set_xticklabels([])
    axes[0, 1].set_xticklabels([])
    axes[0, 1].set_yticklabels([])
    axes[1, 1].set_yticklabels([])
    sb.despine()
    plt.savefig("outputs/c_species_environments.pdf")
    plt.close()
