
# import matplotlib.pyplot as plt

# def plot_per_base_quality(data):
#     """
#     data: list of dicts with keys 'Base' and 'Mean'
#     """
#     x = [d["Base"] for d in data]
#     y = [d["Mean"] for d in data]
#     plt.figure()
#     plt.plot(x, y, marker="o")
#     plt.title("Per-base Sequence Quality (Mean)")
#     plt.xlabel("Base")
#     plt.ylabel("Mean Quality")
#     plt.xticks(rotation=45)
#     plt.tight_layout()
#     return plt.gcf()

# def plot_gc_content(data):
#     x = [d["GC"] for d in data]
#     y = [d["Count"] for d in data]
#     plt.figure()
#     plt.plot(x, y)
#     plt.title("Per-sequence GC Content")
#     plt.xlabel("GC (%)")
#     plt.ylabel("Count")
#     plt.tight_layout()
#     return plt.gcf()

# def plot_length_distribution(data):
#     x = [d["Length"] for d in data]
#     y = [d["Count"] for d in data]
#     plt.figure()
#     plt.plot(x, y, marker="o")
#     plt.title("Sequence Length Distribution")
#     plt.xlabel("Length")
#     plt.ylabel("Count")
#     plt.tight_layout()
#     return plt.gcf()
