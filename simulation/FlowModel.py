import numpy as np
import networkx as nx
import pandas as pd

from simulation import Chain, MarkovChain


class FlowModel:
    def __init__(self, filename, patient_params, simulation_params):
        """
        Constructor for CompartmentNetwork
        total_vol : (L)
        total_flow : (L/s)
        dt : time step in seconds
        """
        self.total_volume = patient_params['TBV']
        self.total_flow = patient_params['CO']
        self.sample_size = simulation_params['sample_size']
        self.nr_steps = simulation_params['nr_steps']
        self.dt = simulation_params['dt']
        self.weibull_shape = simulation_params['weibull_shape']

        self.df = self._read_excel_file(filename, sheetname=patient_params['sheet_name'])
        self.size = self.df.index.name
        self.names = list(self.df.index.values[:self.size])
        # these are given as percentages
        self.flows = np.array(self.df.flow_sum[:self.size].values, dtype=np.float64) / 100 * self.total_flow
        self.volumes = np.array(self.df.volume[:self.size].values, dtype=np.float64) / 100 * self.total_volume
        self.cum_volume = np.cumsum(self.volumes) / np.sum(self.volumes)
        self.particle_volume = self.total_volume / self.sample_size

        # get the rate matrix:
        self.k_matrix = self._get_rate_matrix()
        # convert k_matrix to a graph:
        self._rate_matrix_to_graph()
        self._get_mtts()
        # initialize probability matrix
        self.prob = None

        # initialize chain:
        self.chain = None
        print('Compartmental simulation initialized.')

    def _read_excel_file(self, filename, sheetname):
        df = pd.read_excel(filename, sheet_name=sheetname, engine="openpyxl", index_col=0)
        df.fillna(0, inplace=True)
        return df

    def _get_rate_matrix(self):
        k_matrix = np.array(self.df.values[:self.size, :self.size], dtype=np.float64) / 100 * self.total_flow
        k_matrix /= self.volumes[:, None]
        return k_matrix

    def _get_mtts(self):
        # sort (graph does not have an order necessarily).
        idx = np.array([list(self.G.nodes.keys()).index(name) for name in self.names], dtype=int)
        self.mtt = 1 / np.sum(self.k_matrix, axis=1)[idx]

    def _rate_matrix_to_graph(self):
        # condense flow kinetics in a directed graph with rates k
        adjacency = np.array(self.k_matrix, dtype=[('k', np.float32)])
        self.G = nx.from_numpy_array(adjacency, create_using=nx.DiGraph())
        node_attr = {i: volume for i, volume in enumerate(self.volumes)}
        nx.set_node_attributes(self.G, node_attr, 'V')
        mapping = {i: name for i, name in enumerate(self.names)}
        nx.relabel_nodes(self.G, mapping, copy=False)

    def _graph_to_rate_matrix(self):
        self.k_matrix = nx.to_numpy_array(self.G, weight='k')

    def _get_transition_matrix(self):
        # convert the graph into a transition matrix.
        self._graph_to_rate_matrix()
        self.prob = self.k_matrix * self.dt
        # sort (graph does not have an order necessarily).
        idx = np.array([list(self.G.nodes.keys()).index(name) for name in self.names], dtype=int)
        self.prob = self.prob[idx][:, idx]
        assert(np.array(np.sum(self.prob, axis=1) < np.ones(self.prob.shape[0])).all()), \
            'time step size is too large; leaving probabilities > 1 encountered.'
        # probability of staying
        np.fill_diagonal(self.prob, 1.0 - np.sum(self.prob, axis=1))

    def construct_weibull(self):
        self._get_transition_matrix()
        # Construct jumping process using Weibull distribution
        self.chain = Chain(self.names, self.prob, self.mtt, dt=self.dt, k=self.weibull_shape)

    def construct_markov(self):
        self._get_transition_matrix()
        # Construct markov chain
        self.chain = MarkovChain(self.names, self.prob)


class ExpandFlowModel(FlowModel):
    """
    Inherits from FlowModel. The idea is that we can dynamically add a tumor box.
    In this implementation, the tumor box get added in parallel to the tumor-site from which it 'steals' simulation.
    How much? That is given by the volume fraction (size)
    and the relative simulation density and perfusion of tumor vs tumor-site.
    """
    def __init__(self, filename, patient_params, simulation_params):
        # inherit from base class:
        FlowModel.__init__(self, filename, patient_params, simulation_params)

    def _add_box(self, name, site, blood_volume_fraction):
        idx = self.names.index(site)
        self.size += 1
        self.names.insert(idx + 1, name)

        # adjust volume of original site and added box:
        orig_volume = self.volumes[idx]
        self.volumes[idx] = (1 - blood_volume_fraction) * orig_volume
        self.volumes = np.insert(self.volumes, idx + 1, blood_volume_fraction * orig_volume)
        self.cum_volume = np.cumsum(self.volumes) / np.sum(self.volumes)
        self.G.nodes[site]['V'] = self.volumes[idx]
        self.G.add_node(name, V=self.volumes[idx + 1])
        return idx

    def split_box_parallel(self, name, box_dict):
        site = box_dict['tumor_site']
        volume_fraction = box_dict['tumor_volume_fraction']
        relative_blood_density = box_dict['relative_blood_density']
        relative_perfusion = box_dict['relative_perfusion']

        blood_volume_fraction = volume_fraction * relative_blood_density
        blood_flow_fraction = volume_fraction * relative_perfusion

        assert (blood_volume_fraction < 1.0), \
            'Cannot steal more than 100% of the original simulation volume.'
        assert(blood_flow_fraction < 1.0), \
            'Cannot steal more than 100% of the original flow.'
        idx = self._add_box(name, site, blood_volume_fraction)

        # adjust flow original site and added box:
        orig_flow = self.flows[idx]
        self.flows[idx] = (1 - blood_flow_fraction) * orig_flow
        self.flows = np.insert(self.flows, idx + 1, blood_flow_fraction * orig_flow)

        # adjust network:
        for prev_comp in self.G.predecessors(site):
            orig_rate = self.G.edges[(prev_comp, site)]['k']
            # rate changes as flow since the simulation volume of the predecessor of the site remains equal.
            site_rate = (1 - blood_flow_fraction) * orig_rate
            box_rate = blood_flow_fraction * orig_rate
            self.G.edges[(prev_comp, site)]['k'] = site_rate
            self.G.add_edge(prev_comp, name, k=box_rate)
        for next_comp in self.G.successors(site):
            orig_rate = self.G.edges[(site, next_comp)]['k']
            # Now both the flow and volume are different...
            site_rate = (1 - blood_flow_fraction) / (1 - blood_volume_fraction) * orig_rate
            box_rate = blood_flow_fraction / blood_volume_fraction * orig_rate
            self.G.edges[(site, next_comp)]['k'] = site_rate
            self.G.add_edge(name, next_comp, k=box_rate)

        # the rates have changed so we have to update the rate_matrix and MTTs:
        self._graph_to_rate_matrix()
        self._get_mtts()
