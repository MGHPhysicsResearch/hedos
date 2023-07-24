import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plot_dose_distribution(blood_dose_total, dose_contributions, mean_blood_dose=None):
    x_max1 = 2 * np.percentile(blood_dose_total.dose, 90)
    x_max2 = max([2 * np.percentile(dose, 90) for dose in dose_contributions.values()])
    bins1 = np.linspace(0, x_max1, 100)
    bins2 = np.linspace(0, x_max2, 100)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.hist(blood_dose_total.dose, bins=bins1)
    ax1.axvline(np.mean(blood_dose_total.dose), ymax=1, c='k', linestyle='-',
                label='Simulated mean dose - {:.3f} Gy'.format(np.mean(blood_dose_total.dose)))
    if mean_blood_dose is not None:
        ax1.axvline(mean_blood_dose, ymax=1, c='red', linestyle='--',
                    label='Expected mean dose - {:.3f} Gy'.format(mean_blood_dose))
    for organ, dose in dose_contributions.items():
        ax2.hist(dose, bins=bins2[1:], histtype='step', linewidth=1.5, alpha=0.7, label=organ)

    ax1.set_xlabel('Dose (Gy)')
    ax2.set_xlabel('Dose (Gy)')
    ax1.set_title('Blood dose histogram')
    ax2.set_title('Blood dose contributions')
    ax1.legend()
    ax2.legend()

    plt.show()


def plot_volumes(volume_ref, volume, plot_slice=None, cmap_ref='Greys_r', cmap='Greys_r', scrollable=False):
    """
    Plotting method to visualize 3D volumes.
    Either visualize a slice, or scroll through the entire volume in interactive mode.
    """
    if plot_slice is None:
        plot_slice = np.argmax(np.sum(volume, axis=(0, 1)))

    if scrollable:
        backend = matplotlib.get_backend()
        # you need interactive mode for scrolling, on Mac this works:
        matplotlib.use("QtAgg")
        fig, ax = plt.subplots(1, 1)
        tracker = IndexTracker(ax, volume_ref, volume, plot_slice, cmap_ref, cmap)
        fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
        plt.show()
        # return to original backend.
        matplotlib.use(backend)
    else:
        fig, ax = plt.subplots(1, 1)
        img = ax.imshow(volume_ref[:, :, plot_slice], cmap=cmap_ref)
        c_bar = fig.colorbar(img)
        c_bar.set_label('Treatment dose at slice {} (Gy)'.format(plot_slice))
        ax.imshow(volume[:, :, plot_slice], cmap=cmap, alpha=0.75)
        plt.show()


class IndexTracker(object):
    def __init__(self, ax, X, Y, plot_slice, cmap_ref='Greys', cmap='Greys_r'):
        self.ax = ax
        self.X = X
        self.Y = Y
        self.plot_slice = plot_slice
        _, _, self.slices = X.shape

        self.im1 = ax.imshow(self.X[:, :, self.plot_slice], cmap=cmap_ref)
        self.im2 = ax.imshow(self.Y[:, :, self.plot_slice], cmap=cmap, alpha=0.75)

        c_bar = plt.colorbar(self.im1)
        c_bar.set_label('Treatment dose (Gy)')

        self.update()

    def onscroll(self, event):
        if event.button == 'up':
            self.plot_slice = (self.plot_slice + 1) % self.slices
        else:
            self.plot_slice = (self.plot_slice - 1) % self.slices
        self.update()

    def update(self):
        im1_data = self.im1.to_rgba(self.X[:, :, self.plot_slice], alpha=self.im1.get_alpha())
        im2_data = self.im2.to_rgba(self.Y[:, :, self.plot_slice], alpha=self.im2.get_alpha())

        self.im1.set_data(im1_data)
        self.im2.set_data(im2_data)

        self.ax.set_ylabel('slice %s' % self.plot_slice)
        self.im1.axes.figure.canvas.draw()
        self.im2.axes.figure.canvas.draw()
