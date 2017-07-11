# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: particlestat_tests.py
# Created by: gemusia
# Creation date: 08-07-2017
# Last modified: 11-07-2017 12:45:26
# Purpose:test class for particle statistics (particlestat.py) 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import particlestat as ps
import unittest


#test file: 200_particles
file_path = '/home/gemusia/kody/pp_TCF/MODULES/test_data/'
tf = np.transpose(np.loadtxt(file_path + "200_particles"))

test = ps.Particles(tf[2],tf[0],tf[1],tf[5],tf[3],tf[4],Ux=tf[6])
test1 = ps.Particles([1,2,3],[-1,0,1],[1,2,3],[0.5,0.5,0.5],[1.5,1.5,1.5],[2.5,2.5,2.5])



class ParticleInitTest(unittest.TestCase):

    def test_input0(self):
        try:
            ps.Particles(tf[2],tf[0],tf[1],tf[5],tf[3],tf[4])
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_input1(self):
        with self.assertRaises(TypeError):
            ps.Particles(tf[2],tf[0],tf[1],tf[5],tf[3],tf[4],tf[6])

    def test_input2(self):
        try:
            ps.Particles(tf[2],tf[0],tf[1],tf[5],tf[3],tf[4],Ux=tf[6])
        except ValueError:
            self.fail("an unexpected ValueError")


class ParticleStatTest(unittest.TestCase):

    def test_mean0(self):
        with self.assertRaises(TypeError):
            test.pmean()

    def test_mean1(self):
        try:
            test.pmean("x")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_mean2(self):
        try:
            test.pmean("Ux")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_mean3(self):
        self.assertTrue(np.allclose(test.pmean("y")[0][0:16],test.pmean("y")[1][0:16],atol=0.2))

    def test_pvar(self):
        try:
            test.pvar("Ux")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pstd(self):
        try:
            test.pvar("Ux")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pcov(self):
        try:
            test.pcov("Vx","Vy")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pke(self):
        try:
            test.pke("Vx","Vy","Vz")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_mean_symm(self):
        try:
            test.stat_symm("pmean","symm","Ux")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pvar_symm(self):
        try:
            test.stat_symm("pvar","symm","Ux")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pstd_symm(self):
        try:
            test.stat_symm("pstd","symm","Ux")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pcov_symm(self):
        try:
            test.stat_symm("pcov","symm","Vx","Vy")
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_pke_symm(self):
        try:
            test.stat_symm("pke","symm","Vx","Vy","Vz")
        except ValueError:
            self.fail("an unexpected ValueError")
#******************************************
#call of the program
#******************************************

if __name__=='__main__':
    unittest.main()

