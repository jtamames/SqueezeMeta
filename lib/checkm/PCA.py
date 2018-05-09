###############################################################################
#
# PCA.py - Principal component analysis
#
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

    def pcaFile(self, dataMatrixFile, fraction=0.90, bCenter=True, bScale=True):
        """ Perform PCA on TSV file with row and column headers.
              Rows are data points, columns are variables."""
        data = np.array([])
        names = np.array([])
        bHeader = True
        rows = 0
        for line in open(dataMatrixFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')

            rows += 1
            cols = len(lineSplit) - 1

            names = np.append(names, [lineSplit[0]])
            data = np.append(data, [float(x) for x in line.split('\t')[1:]])

        data = np.reshape(data, (rows, cols))

        # do the PCA and extract the scores
        self.pcaMatrix(data, fraction, bCenter, bScale)

        return names, self.pc(), self.variance

    def pcaMatrix(self, A, fraction=0.90, bCenter=True, bScale=True):
        """Perform PCA of a data matrix. Rows are data points, columns are variables."""
        assert 0 <= fraction <= 1

        if bCenter:
            Center(A, bScale)

        self.U, self.d, self.Vt = np.linalg.svd(A, full_matrices=False)
        assert np.all(self.d[:-1] >= self.d[1:])  # sorted

        self.variance = self.d ** 2
        self.variance /= np.sum(self.variance)

        self.sumvariance = np.cumsum(self.variance)

        self.npc = np.searchsorted(self.sumvariance, fraction) + 1
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
    def __init__(self, A, axis=0, scale=True):
        self.mean = A.mean(axis=axis)
        A -= self.mean
        if scale:
            std = A.std(axis=axis)
            self.std = np.where(std, std, 1.)
            A /= self.std
        else:
            self.std = np.ones(A.shape[-1])
        self.A = A

    def uncenter(self, x):
        return np.dot(self.A, x * self.std) + np.dot(x, self.mean)
