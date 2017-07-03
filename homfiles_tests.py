# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat_tests.py
# Created by: gemusia
# Creation date: 23-06-2017
# Last modified: 03-07-2017 16:26:37
# Purpose: test of module homfiles.py
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import homstat as hs
import homfiles as hf
import unittest


test_set = hf.ChannelFields(128,129,128,2501,2505)
test_set1 = hf.ChannelFields(128,129,128,2501,2502)

#******************************************
#  __init__  tests
#******************************************
class ChannelTest(unittest.TestCase):

    def test_input(self):
        self.assertSetEqual(set(test_set.__dict__.keys()),
                set(['nFiles','fStart','fEnd','fNames','field_2501','field_2502','field_2503','field_2504','field_2505']))



#******************************************
# test for statistics computations
#******************************************

#mean
class ChannelTest_mean(unittest.TestCase):


    def test_mean_symm(self):
        thmean = test_set1.mean_symmT()
        print thmean[1]
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thmean[i],np.zeros(65)))




'''

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
        thstd = tst_std.hstd()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.full(3,0.5)))


    def test_std_symm0(self):
        thstd = test_128.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(64)))

    def test_std_symm1(self):
        thstd = test_64_32.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(16)))

    def test_std_symm2(self):
        thstd = test_128ones.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(64)))


    def test_std_symm3(self):
        thstd = test_32ones.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.zeros(16)))


    def test_std_symm4(self):
        thstd = tst_std.hstd_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.full(2,0.5)))





#cor - standard deviation
class ChannelTest_cor(unittest.TestCase):


    def test_cor0(self):
        for i in range(3):
            for j in range(3):
                thcor = test_128.hcor(U_dict[i],U_dict[j])
                self.assertTrue(np.allclose(thcor[1],np.zeros(128)))

    def test_cor1(self):
        for i in range(3):
            thstd = tst_cor.hstd()
            for j in range(3):
                thcor = tst_cor.hcor(U_dict[i],U_dict[j])
                if i<>j:
                    cor_coeff= np.divide(thcor[1],thstd[i+1]*thstd[j+1])
                    self.assertTrue(np.allclose(cor_coeff,np.ones(3)),"i = %d, j = %d" %(i,j))
                else:
                    self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])))

    def test_cor2(self):
        for i in range(3):
            thstd = tst1_cor.hstd()
            for j in range(3):
                thcor = tst1_cor.hcor(U_dict[i],U_dict[j])
                cor_coeff= np.divide(thcor[1],thstd[i+1]*thstd[j+1])
                if (i,j) in {(0,1),(1,0),(0,2),(2,0)}:
                    self.assertTrue(np.allclose(cor_coeff,np.full(3,-1.0)),"i = %d, j = %d" %(i,j))
                elif (i,j) in {(1,2),(2,1)}:
                    self.assertTrue(np.allclose(cor_coeff,np.full(3,1.0)),"i = %d, j = %d" %(i,j))
                else:
                    self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])))



    def test_cor_symm0(self):
        for i in range(3):
            for j in range(3):
                thcor = test_128.hcor_symm(U_dict[i],U_dict[j])
                self.assertTrue(np.allclose(thcor[1],np.zeros(64)))

    def test_cor_symm1(self):
        for i in range(3):
            thstd = tst_cor.hstd_symm()
            for j in range(3):
                thcor = tst_cor.hcor_symm(U_dict[i],U_dict[j])
                if (i,j) in {(0,1),(1,0),(1,2),(2,1)}:
                    self.assertTrue(np.allclose(thcor[1],np.zeros(2)),"i = %d, j = %d" %(i,j))
                elif (i,j) in {(0,2),(2,0)}:
                    cor_coeff= np.divide(thcor[1],thstd[i+1]*thstd[j+1])
                    self.assertTrue(np.allclose(cor_coeff,np.ones(2)),"i = %d, j = %d" %(i,j))
                else:
                    self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])))

    def test_cor_symm2(self):
        for i in range(3):
            thstd = tst1_cor.hstd_symm()
            for j in range(3):
                thcor = tst1_cor.hcor_symm(U_dict[i],U_dict[j])
                cor_coeff= np.divide(thcor[1],thstd[i+1]*thstd[j+1])
                if (i,j) in {(0,2),(2,0)}:
                    self.assertTrue(np.allclose(cor_coeff,np.full(2,-1.0)),"i = %d, j = %d" %(i,j))
                elif (i,j) in {(0,1),(1,0),(1,2),(2,1)}:
                    self.assertTrue(np.allclose(cor_coeff,np.zeros(2)),"i = %d, j = %d" %(i,j))
                else:
                    self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])))

    #test of cerreation coefficient computation on real data
    # correlation (i,i) is compared to square of standard deviation
    def test_cor_data(self):
	for k in Data_dict.keys():
	    for i in range(3):
		thstd = Data_dict[k].hstd()
		thcor = Data_dict[k].hcor(U_dict[i],U_dict[i])
		self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])))


    def test_cor_data_symm(self):
	for k in Data_dict.keys():
	    for i in range(3):
		thstd = Data_dict[k].hstd_symm()
		thcor = Data_dict[k].hcor_symm(U_dict[i],U_dict[i])
		self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1]),atol=1e-04))

    #symmetrisation of cross-correlations is roughly tested
    def test_cor_data_symm(self):
	for k in Data_dict.keys():
	    for i in range(3):
                for j in range(3):
		    thcor = Data_dict[k].hcor(U_dict[i],U_dict[j])
		    thcor_symm = Data_dict[k].hcor_symm(U_dict[i],U_dict[j])
		    self.assertTrue(np.allclose(thcor[1][len(thcor_symm[1])-1:],np.flipud(thcor_symm[1]),atol=0.01))


#TODO - write tests for kinetic energy



'''

#******************************************
#call of the program
#******************************************

if __name__=='__main__':
    unittest.main()

