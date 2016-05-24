import numpy as np
import GPy
import csv
import matplotlib

code_dir = os.path.expanduser("~/code/ebola_lofa/")
data_dir = os.path.expanduser("~/Data/Ebola/Lofa/")

admissions = np.genfromtxt(code_dir + "/data/admissions.csv", delimiter = ',', skip_header = 1)
deaths = np.genfromtxt(code_dir + "/data/deaths.csv", delimiter = ',', skip_header = 1)

kernel = GPy.kern.RBF(input_dim = 1, variance = 1., lengthscale = 1.)
m_adm = GPy.models.GPRegression(np.transpose(np.atleast_2d(admissions[:,0])),
                                np.transpose(np.atleast_2d(admissions[:,1])),
                                kernel)
m_adm.Gaussian_noise.variance.constrain_bounded(0,2)

m_adm.optimize()
print(m_adm)

fig_adm = m_adm.plot()
fig_adm.figure.savefig(data_dir + "gp_admissions.pdf")

m_dea = GPy.models.GPRegression(np.transpose(np.atleast_2d(deaths[:,0])),
                            np.transpose(np.atleast_2d(deaths[:,1])),
                            kernel)

m_dea.Gaussian_noise.variance.constrain_bounded(0,2)

m_dea.optimize()
print(m_dea)

fig_dea = m_dea.plot()
fig_dea.figure.savefig(data_dir + "gp_deaths.pdf")

