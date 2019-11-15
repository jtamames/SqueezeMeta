###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import numpy as np


class PCA:
    """http://stackoverflow.com/questions/1730600/principal-component-analysis-in-python"""
    def __init__(self):
        pass

    def pca_matrix(self, A, num_components, bCenter=True, bScale=True):
        """Perform PCA of a data matrix. Rows are data points, columns are variables."""

        if bCenter:
            Center(A, bScale=bScale)

        self.U, self.d, self.Vt = np.linalg.svd(A, full_matrices=False)
        assert np.all(self.d[:-1] >= self.d[1:])  # sorted

        self.variance = self.d ** 2
        self.variance /= np.sum(self.variance)
        self.sumvariance = np.cumsum(self.variance)

        self.npc = min(num_components, len(self.variance))
        self.dinv = np.array([1 / d if d > self.d[0] * 1e-6  else 0
                                for d in self.d])

        return self.pc(), self.variance

    def pc(self):
        """ e.g. 1000 x 2 U[:, :npc] * d[:npc], to plot etc. """
        n = self.npc
        return self.U[:, :n] * self.d[:n]

    def vars_pc(self, x):
        n = self.npc
        return self.d[:n] * np.dot(self.Vt[:n], x.T).T  # 20 vars -> 2 principal

    def pc_vars(self, p):
        n = self.npc
        return np.dot(self.Vt[:n].T, (self.dinv[:n] * p).T) .T  # 2 PC -> 20 vars

    def pc_obs(self, p):
        n = self.npc
        return np.dot(self.U[:, :n], p.T)  # 2 principal -> 1000 obs

    def obs_pc(self, obs):
        n = self.npc
        return np.dot(self.U[:, :n].T, obs) .T  # 1000 obs -> 2 principal

    def obs(self, x):
        return self.pc_obs(self.vars_pc(x))  # 20 vars -> 2 principal -> 1000 obs

    def vars(self, obs):
        return self.pc_vars(self.obs_pc(obs))  # 1000 obs -> 2 principal -> 20 vars


class Center:
    """http://stackoverflow.com/questions/1730600/principal-component-analysis-in-python"""
    """ A -= A.mean() /= A.std(), inplace -- use A.copy() if need be
        uncenter(x) == original A . x
    """
    def __init__(self, A, axis=0, bScale=True):

        self.mean = A.mean(axis=axis)
        A -= self.mean

        if bScale:
            std = A.std(axis=axis)
            self.std = np.where(std, std, 1.)
            A /= self.std
        else:
            self.std = np.ones(A.shape[-1])

        self.A = A

    def uncenter(self, x):
        return np.dot(self.A, x * self.std) + np.dot(x, self.mean)
