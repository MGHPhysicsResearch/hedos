import numpy as np
import matplotlib.pyplot as plt
import time


class TemporalDistribution:
    def __init__(self, model):
        self.model = model
        self.ttd = {}
        self.rtd = {}
        self.path = None
        self.tv = None

    def generate_from_markov(self):
        """
        Generate a temporal distribution from a pre-built markov chain.
        """
        t = time.process_time()
        compartment_id = self.model.cum_volume.searchsorted(
            np.random.uniform(size=self.model.sample_size)).astype(np.uint8)
        self.path = self.model.chain.walk(self.model.nr_steps, compartment_id)
        print(f'Time to generate simulation distribution: {time.process_time()-t:.6f} seconds')

    def generate_from_weibull(self):
        """
        Generate a temporal distribution from a pre-built chain using Weibull distribution.
        """
        start_time = time.process_time()
        compartment_id = self.model.cum_volume.searchsorted(
            np.random.uniform(size=self.model.sample_size)).astype(np.uint8)
        self.path = self.model.chain.walk_v1(self.model.nr_steps, compartment_id)
        # self.path = self.model.chain.walk_v2(self.model.nr_steps, compartment_id)
        print(f'Time to generate temporal distribution: {time.process_time()-start_time:.2f} seconds')

    def temporal_volume(self):
        """
        Calculate volume changes in time. # of BP x # of time-steps
        """
        start_time = time.process_time()
        self.tv = np.apply_along_axis(lambda x: np.bincount(x, minlength=256), axis=0, arr=self.path)
        self.tv = self.tv[:len(self.model.names)] / self.model.sample_size
        print(f'Time to get temporal volumes: {time.process_time() - start_time:.2f} seconds')

    def save(self, f_name):
        """
        Save simulation path for potential re-use
        """
        np.save(f_name, self.path)

    def load(self, f_name):
        """
        Load simulation distribution
        """
        self.path = np.load(f_name)

    def _particle_entry_exit(self, compartment_id):
        # which rows (=particle ids) pass given compartment at some point during simulated time:
        rows = np.amax(self.path == compartment_id, axis=1)
        # for these particle ids, when is it in the compartment:
        in_comp = np.array(self.path[rows] == compartment_id, dtype=np.int32)
        # when does it enter and leave the compartment?:
        diff = np.concatenate([in_comp[:, 0][:, None], np.diff(in_comp, axis=1), -in_comp[:, -1][:, None]], axis=1)
        particle_id, t_entry = np.where(diff == 1)
        t_exit = np.where(diff == -1)[1]
        return in_comp, particle_id, t_entry, t_exit

    def _get_time_distributions(self, name, nr_particles_passed):
        """
        Due to finite width of time window, this is biased toward shorter transition/recurrence times as their
        probability of falling within the window in greater. I think the correct term for this is "right censoring".
        To quantify this a bit, we print the percentage of transit/recurrence times that were contained
        in the simulated time window.
        For this we need an estimate of particles that have past the compartment during the simulation.
        This is given by the "nr_particles_passed" parameter.
        """
        # find the time indices where the particle enters and exits the compartment.
        # get transition and recurrence times by subtraction.
        compartment_id = self.model.names.index(name)
        in_comp, particle_id, t_entry, t_exit = self._particle_entry_exit(compartment_id)
        # Calculate tts
        ttd = (t_exit - t_entry) * self.model.dt
        # beginning and tail are potentially cut off due to specific time window, discard these.
        self.ttd[name] = ttd[in_comp[particle_id, 0] + in_comp[particle_id, -1] == 0]
        # Calculate rts, deleting entries that do not correspond to the same particle.
        rtd = (t_entry[1:] - t_exit[:-1]) * self.model.dt
        self.rtd[name] = rtd[np.where(particle_id[1:] == particle_id[:-1])]

        print('Only the shortest {:.3f}% of transit times captured in {:}, MTT therefore underestimation.'.format(
            self.ttd[name].size / nr_particles_passed * 100, name))
        print('Only the shortest {:.3f}% of recurrence times captured in {:}, MRT therefore underestimation.'.format(
            self.rtd[name].size / nr_particles_passed * 100, name))

    def _transition_recurrence_time(self, expected_nr_particles, names=None):
        """
        Calculate transition time tt-distribution (ttd), and recurrence time rt-distribution (rtd)
        """
        if names is None:
            names = self.model.names

        start_time = time.process_time()
        for name in names:
            nr_particles_passed = expected_nr_particles[self.model.names.index(name)]
            self._get_time_distributions(name, nr_particles_passed)
        print(f'Time to get transition times: {time.process_time() - start_time:.6f} seconds')

    # these are just some plotting functions to get a feel for the generated spatiotemporal distribution...
    def _plot_hist(self, names, time_distribution, name_of_mean):
        _, axes = plt.subplots(nrows=len(names), ncols=1, figsize=(6, 4 * len(names)))
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        for ax, name in zip(axes, names):
            v_max = 3 * np.percentile(time_distribution[name], 90)
            n_bins = int(v_max // self.model.dt)
            ax.hist(time_distribution[name], bins=np.linspace(0, v_max, n_bins),
                    density=True, label=name)
            ax.legend()
            ax.set_title(name_of_mean + ' = {:.3f}'.format(np.mean(time_distribution[name])))
        plt.xlabel('time (s)')
        plt.show()

    def plot_time_distributions(self, names):
        """
        This plots both the simulated transit time distribution and the recurrence time distribution.
        Note that both are right-censored; their mean will therefore be an underestimation.
        """
        assert(isinstance(names, list)), '"names" should be a list.'
        # obtain the expected nr of particles crossing through each organ in the given time window.
        volume_passed = self.model.flows * self.path.shape[1] * self.model.dt
        expected_nr_particles = volume_passed / self.model.particle_volume
        # calculate transit and recurrence time distributions.
        self._transition_recurrence_time(expected_nr_particles, names)

        self._plot_hist(names, time_distribution=self.ttd, name_of_mean='MTT')
        self._plot_hist(names, time_distribution=self.rtd, name_of_mean='MRT')

    def plot_inflow_outflow(self, names):
        """
        This plots both the flow into and out off a compartment.
        Clearly in equilibrium these should be the same and equal the intended compartmental flow.
        """
        for name in names:
            compartment_id = self.model.names.index(name)
            in_comp, particle_id, t_entry, t_exit = self._particle_entry_exit(compartment_id)
            ti, count = np.unique(t_entry, return_counts=True)
            plt.plot(ti[1:-1] * self.model.dt, count[1:-1] * self.model.particle_volume * 1000 / self.model.dt,
                     label=name + ' -- inflow')
            ti, count = np.unique(t_exit, return_counts=True)
            plt.plot(ti[1:-1] * self.model.dt, count[1:-1] * self.model.particle_volume * 1000 / self.model.dt,
                     label=name + ' -- outflow')
        plt.xlim([0, self.path.shape[1] * self.model.dt])
        plt.legend()
        plt.xlabel('time (s)')
        plt.ylabel('flow (mL/s)')
        plt.show()

    def plot_volumes_over_time(self):
        """
        This plots the (simulation) volumes of the all compartments over time.
        Again, in equilibrium, these lines should be straight (with some stochastic noise)
        and reflect the intended simulation volume in each compartment.
        """
        if self.tv is None:
            self.temporal_volume()
        ti = np.arange(self.tv.shape[1]) * self.model.dt
        for i, name in enumerate(self.model.names):
            plt.plot(ti, self.tv[i], label=name)
        plt.legend()
        plt.xlabel('time (s)')
        plt.show()

    def plot_normalized_volume_over_time(self, names):
        """
        This zooms in on one or multiple volumes and plots them normalized.
        Hence, in equilibrium, these lines should wiggle about the value 1.
        """
        if self.tv is None:
            self.temporal_volume()
        # normalize volumes:
        normalized_volumes = self.model.volumes / np.sum(self.model.volumes)
        print('Sum target volumes = {:.3f}'.format(np.sum(normalized_volumes)))
        final_volumes = self.tv[:, -1]
        print('Sum final volumes = {:.3f}'.format(np.sum(final_volumes)))

        _, axes = plt.subplots(nrows=len(names), ncols=1, figsize=(6, 4 * len(names)))
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        for ax, name in zip(axes, names):
            idx = self.model.names.index(name)
            ax.plot(np.arange(self.tv[idx].size) * self.model.dt, self.tv[idx] / normalized_volumes[idx], label=name)
            ax.legend()
            print('{} -- Fraction of target volume = {:.3f}'.format(name, final_volumes[idx] / normalized_volumes[idx]))
        plt.xlabel('time (s)')
        plt.suptitle('Normalized volume over time')
        plt.show()

    def plot_final_blood_volumes(self):
        """
        This plots the final volumes + reference values.
        Perhaps better to take an average of the final x steps to get rid of some stochasticity.
        """
        bins = self.model.size
        fig, ax = plt.subplots(1, 1, figsize=(5, 8), dpi=300)
        ax.hist(self.path[:, -1], bins=bins, range=[-0.5, bins - 0.5],
                weights=[1/self.model.sample_size] * self.model.sample_size, density=False,
                label='simulation', orientation='horizontal', rwidth=0.7)
        ax.plot(self.model.volumes / np.sum(self.model.volumes), np.arange(bins), 'r*', label='reference')
        ax.set_yticks(list(range(0, bins)), self.model.names)
        ax.set_xticks([0, 0.05, 0.1, 0.15], [0, 5, 10, 15])
        ax.set_xlabel('Percentage of total simulation volume')
        ax.invert_yaxis()
        plt.subplots_adjust(bottom=0.4)
        plt.legend()
        plt.tight_layout()
        plt.show()


