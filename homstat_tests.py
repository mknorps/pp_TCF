# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat_tests.py
# Created by: gemusia
# Creation date: 23-06-2017
# Last modified: 29-06-2017 15:37:49
# Purpose: test of module homstat.py
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import homstat as hs
import unittest

#******************************************
#useful input examples and constants
#******************************************

#real DNS dimensions
zeros128 = np.zeros(shape=(128,128,128))
test_128 = hs.Channel(zeros128,zeros128,zeros128,128,128,128)

ones128 = np.ones(shape=(128,128,128))
test_128ones = hs.Channel(ones128,ones128,ones128,128,128,128)


#real LES dimensions
zeros64_32 = np.zeros(shape=(32,32,64))
test_64_32 = hs.Channel(zeros64_32,zeros64_32,zeros64_32,32,32,64)

ones_32 = np.ones(shape=(32,32,64))
test_32ones = hs.Channel(ones_32,ones_32,ones_32,32,32,64)

#nodes of Chebyshev
y_NP64 = np.loadtxt("y_nodes.txt")



#******************************************
#input tests
#******************************************
class ChannelTest(unittest.TestCase):

    def test_input0(self):
        with self.assertRaises(ValueError):
            hs.Channel([],[],[],0,0,0)

    def test_input1(self):
        try:
            hs.Channel([[[1]]],[[[2]]],[[[3]]],1,1,1)
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_input2(self):
        with self.assertRaises(ValueError):
            hs.Channel([[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],1,2,1)

    def test_input3(self):
        try:
            hs.Channel([[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],2,2,1)
        except ValueError:
            self.fail("an unexpected ValueError")

    def test_inputi4(self):
        with self.assertRaises(ValueError):
            hs.Channel([[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],2,1,1)

    def test_input5(self):
        try:
            test_128
        except ValueError:
            self.fail("an unexpected ValueError")




#wall normal direction node computation tests
class ChannelTest_nondim(unittest.TestCase):

    def test_ynondim(self):
        self.assertTrue(np.allclose(test_128.y_nondim(),y_NP64))





#******************************************
# test for statistics computations
#******************************************

#mean
class ChannelTest_mean(unittest.TestCase):

    def test_mean0(self):
        thmean = test_128.hmean()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thmean[i],np.zeros(128)))


    def test_mean1(self):
        thmean = test_64_32.hmean()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thmean[i],np.zeros(32)))


    def test_mean_symm0(self):
        thmean = test_128.hmean_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thmean[i],np.zeros(65)))


    def test_mean_symm1(self):
        thmean = test_64_32.hmean_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thmean[i],np.zeros(17)))

    def test_mean_symm2(self):
        thmean = test_128ones.hmean_symm()
        self.assertTrue(np.allclose(thmean[2],np.zeros(65)))
        for i in {1,3}:
            self.assertTrue(np.allclose(thmean[i],np.ones(65)))


    def test_mean_symm3(self):
        thmean = test_32ones.hmean_symm()
        self.assertTrue(np.allclose(thmean[2],np.zeros(17)))
        for i in {1,3}:
            self.assertTrue(np.allclose(thmean[i],np.ones(17)))


#std - standard deviation
class ChannelTest_std(unittest.TestCase):

    def test_std0(self):
        thstd = test_128.hstd()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(128)))

    def test_std1(self):
        thstd = test_64_32.hstd()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(32)))

    def test_std2(self):
        tst = [[[1,0,0],[1,0,0],[1,0,0]],[[1,0,0],[1,0,0],[1,0,0]],[[1,0,0],[1,0,0],[1,0,0]]]
        tst_channel = hs.Channel(tst,tst,tst,3,3,3)
        thstd = tst_channel.hstd()
        print thstd
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(3)))

    def test_std3(self):
        thstd = test_64_32.hstd()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(32)))

    def test_std_symm0(self):
        thstd = test_128.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(65)))


    def test_std_symm1(self):
        thstd = test_64_32.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(17)))

    def test_std_symm2(self):
        thstd = test_128ones.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(65)))


    def test_std_symm3(self):
        thstd = test_32ones.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(17)))


#******************************************
#tests for functions in module homstat
#******************************************
class HomstatTest_functions(unittest.TestCase):
 
    def test_cheb0(self):
        for x in {32,64,128}:
            self.assertAlmostEqual(hs.ChebZeros(x)(0), 1.0)
  
    def test_cheb1(self):
        for x in {32,64,128}:
            self.assertAlmostEqual(hs.ChebZeros(x)(x), -1.0)

    def test_cheb2(self):
        for x in {32,64,128}:
            self.assertAlmostEqual(hs.ChebZeros(x)(x/2), 0.0)

    def test_cheb3(self):
        for x in {32,64,128}:
            self.assertAlmostEqual(hs.ChebZeros(x)(x/4), np.cos(np.pi/4))


    def test_symm0(self):
        self.assertTrue(np.allclose(hs.symm("symm",np.ones(33)),np.ones(17)))

    def test_symm1(self):
        self.assertTrue(np.allclose(hs.symm("asymm",np.ones(33)),np.zeros(17)))

    def test_symm2(self):
        self.assertTrue(np.allclose(hs.symm("symm",np.full(14,3.14)),np.full(8,3.14)))

    def test_symm3(self):
        self.assertTrue(np.allclose(hs.symm("asymm",np.full(32,3.1415)),np.zeros(17)))





#******************************************
#call of the program
#******************************************

if __name__=='__main__':
    unittest.main()

