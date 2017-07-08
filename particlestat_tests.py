# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: particlestat_tests.py
# Created by: gemusia
# Creation date: 08-07-2017
# Last modified: 08-07-2017 16:02:16
# Purpose:test class for particle statistics (particlestat.py) 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import particlestat as ps
import unittest


#test file: 200_particles
tf = np.transpose(np.loadtxt("200_particles"))

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


class ParticleMeanTest(unittest.TestCase):

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


#******************************************
#call of the program
#******************************************

if __name__=='__main__':
    unittest.main()

