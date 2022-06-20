import networkx as nx
import numpy as np
import pandas as pd
import pydtmc  as dtmc

from blooddvh import TimeDependentMarkovChain


class CompartmentModel:
    __slots__ = ['fin', 'df', 'size', 'name', 'total_volume',
                 'cardiac_output', 'G', 'prob', 'volume', 'markov', 'scales',
                 'markov_weibull']
    def __init__(self, f_name, s_name, vol=5.3, cardiac=6.5, resolution=60):
        """
        Constructor for CompartmentModel.

        Parameters
        ----------
        f_name : str
            The name of the Excel file with the desired blood distribution.
        s_name : str
            The name of the sheet with the desired blood distribution from
            the Excel file, `f_name`.
        vol : float/int, optional
            The total blood volume (L), e.g., ICRP male 5.3 L.
        cardiac : float/int, optional
            The cardiac output (L), e.g., 6.5 L/min.
        resolution : float/int, optional
            The seconds for the cardiac output, 60 means 1 min, i.e.,
            6.5 L/min.

        Returns
        -------
        N/A

        """
        self.fin = {'file': f_name, 'sheet' : s_name}
        self.df = pd.read_excel(self.fin['file'], sheet_name=self.fin['sheet'], engine='openpyxl')
        self.df.fillna(0, inplace=True)
        self.size = int(self.df.columns[0])
        self.name = [c for c in self.df.columns[1:self.size+1]]
        self.volume = np.cumsum(self.df.volume[0:self.size].values) / np.sum(self.df.volume[0:self.size].values)
        self.scales = np.zeros(self.size)
        self.total_volume = vol

        # Flow per sec (L/s)
        self.cardiac_output = cardiac/resolution

        # Create the network
        self.G = nx.DiGraph()

        # 1. Add nodes
        for _,c in enumerate(self.name):
            self.G.add_node(c)

        # 2. Create edges with transition probabilities
        self.prob = np.array(self.df.values[0:self.size, 1:self.size+1], dtype=np.float64)
        for row in range(self.size):
            # Normalize to prevent non-100
            # self.df.flow_sum = sum(self.prob[row])
            # Calculate the probability of leaving
            for col in range(self.size):
                self.prob[row, col] *= self.cardiac_output / 100.0
                self.prob[row, col] /= self.total_volume * self.df.volume[row] / 100.0
            # Determine the probability of staying
            self.prob[row,row] = 1.0 - sum(self.prob[row])
            self.scales[row] = 0.01 * self.total_volume * self.df.volume[row]
            self.scales[row] /= (cardiac/resolution * 0.01 * self.df.flow_sum[row])
            # Create a network edge for non-zero transition
            for col in range(self.size):
                if self.prob[row, col] > 0:
                    self.G.add_edge(self.name[row], self.name[col], weight=self.prob[row,col])

        # 3. Create a Markov chain
        self.markov = dtmc.MarkovChain(self.prob, self.name)

        # 4. Create a time-dependent Markov chain using a Weibull distribution
        shapes_two = 2.0 * np.ones(self.size)
        self.markov_weibull = TimeDependentMarkovChain(self.prob, self.scales, shapes_two)

    def describe(self, options):
        """
        Print specific output.

        Parameters
        ----------
        option : list[str]
            How to determine what to print. Can either be `n_nodes`,
            `n_edges`, `all_nodes`, `all_edges`, `in_degree`, or `out_degree`.

        Returns
        -------
        N/A

        """
        for option in options:
            if option == 'n_nodes':
                print('Total number of nodes: ', int(self.G.number_of_nodes()))
            elif option == 'n_edges':
                print('Total number of edges: ', int(self.G.number_of_edges()))
            elif option == 'all_nodes':
                print('List of all nodes: ', list(self.G.nodes()))
            elif option == 'all_edges':
                print('List of all edges: ', list(self.G.edges()))
            elif option == 'in_degree':
                print('In-degree for all nodes: ', dict(self.G.in_degree()))
            elif option == 'out_degree':
                print('Out degree for all nodes: ', dict(self.G.out_degree))
            else:
                raise Exception('Invalid `option` parameter.')
