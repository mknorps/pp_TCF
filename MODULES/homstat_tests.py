# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat_tests.py
# Created by: gemusia
# Creation date: 23-06-2017
# Last modified: 11-07-2017 12:43:07
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

zeros129 = np.zeros(shape=(128,129,128))
test_129 = hs.Channel(zeros129,zeros129,zeros129,128,129,128)

ones128 = np.ones(shape=(128,128,128))
test_128ones = hs.Channel(ones128,ones128,ones128,128,128,128)


#real LES dimensions
zeros64_32 = np.zeros(shape=(32,32,64))
test_64_32 = hs.Channel(zeros64_32,zeros64_32,zeros64_32,32,32,64)

ones_32 = np.ones(shape=(32,32,64))
test_32ones = hs.Channel(ones_32,ones_32,ones_32,32,32,64)

# input examples for standard deviation and correlation testing
tst = [[[1,0,0,1],[1,0,0,1],[1,0,0,1]],[[1,0,0,1],[1,0,0,1],[1,0,0,1]]]
tst1 = [[[0,1,1,0],[0,1,1,0],[0,1,1,0]],[[0,1,1,0],[0,1,1,0],[0,1,1,0]]]

tst_std = hs.Channel(tst,tst,tst,2,3,4)
tst_cor = hs.Channel(tst,tst,tst,2,3,4)
tst1_cor = hs.Channel(tst1,tst,tst,2,3,4)

U_dict = {0:"Ux",1:"Uy",2:"Uz"}



#******************************************
# real data input; 
#******************************************

# WARNING - data from fortran file was written in different coordinates
# so we have to take Ux as Uy, Uy as Uz and Uz as Ux (permutation! )
file_path = '/home/gemusia/kody/pp_TCF/MODULES/test_data/'
fName  = file_path + "upp_"     #DNS
fNamef = file_path + "uppf_"    #LES (a priori DNS)
upp   = [[]]*3
uppf  = [[]]*3

for s in U_dict.keys():
    with open(fName + U_dict[s] + "_2505",'r') as fupp:
            # input data is in a matrix, and we first have to reshape 
            # it to be a 3D field
            # then directions must be permuted (np.transpose)
	    upp[np.mod(s+1,3)]=np.transpose(np.loadtxt(fupp).reshape(129,128,128),axes=(2,0,1))
    with open(fNamef + U_dict[s] + "_2505",'r') as fuppf:
	    uppf[np.mod(s+1,3)]=np.transpose(np.loadtxt(fuppf).reshape(33,64,32),axes=(2,0,1))

#initialize Channel objects
U_DNS  = hs.Channel(upp[0],upp[1],upp[2],128,129,128) 
U_LES = hs.Channel(uppf[0],uppf[1],uppf[2],32,33,64)


Data_dict = {0:U_DNS,1:U_LES}

#nodes of Chebyshev
y_NP64 = np.loadtxt(file_path + "y_nodes.txt")


#******************************************
#  __init__  tests
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

    def test_ycenters(self):
        self.assertTrue(np.allclose(test_129.y_centers(),y_NP64))

    def test_ynondim1(self):
        self.assertTrue(np.allclose(tst_std.y_nondim(),np.array([0,150])))




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
            self.assertTrue(np.allclose(thmean[i],np.zeros(64)))


    def test_mean_symm1(self):
        thmean = test_64_32.hmean_symm()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thmean[i],np.zeros(16)))

    def test_mean_symm2(self):
        thmean = test_128ones.hmean_symm()
        self.assertTrue(np.allclose(thmean[2],np.zeros(64)))
        und = test_128ones.utau()
        for i in {1,3}:
            self.assertTrue(np.allclose(thmean[i],np.full(64,1.0/und)))


    def test_mean_symm3(self):
        thmean = test_32ones.hmean_symm()
        self.assertTrue(np.allclose(thmean[2],np.zeros(16)))
        und = test_32ones.utau()
        for i in {1,3}:
            self.assertTrue(np.allclose(thmean[i],np.full(16,1.0/und)))




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
        und = tst_std.utau()
        for i in {1,2,3}:
            self.assertTrue(np.allclose(thstd[i],np.full(2,0.5/und)))





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
        und = tst_cor.utau()
        for i in range(3):
            thstd = tst_cor.hstd_symm()
            for j in range(3):
                thcor = tst_cor.hcor_symm(U_dict[i],U_dict[j])
                if (i,j) in {(0,1),(1,0),(1,2),(2,1)}:
                    self.assertTrue(np.allclose(thcor[1],np.zeros(2)),"i = %d, j = %d" %(i,j))
                elif (i,j) in {(0,2),(2,0)}:
                    cor_coeff= np.divide(thcor[1],thstd[i+1]*thstd[j+1])
                    self.assertTrue(np.allclose(cor_coeff,np.full(2,1.0)),"i = %d, j = %d" %(i,j))
                else:
                    self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])),"diagonal component %d" %j)

    def test_cor_symm2(self):
        und = tst_cor.utau()
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
                    self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1])), "diagonal component %d" %j)

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
		self.assertTrue(np.allclose(thcor[1],map(lambda x: x**2,thstd[i+1]),atol=1e-01))

    #symmetrisation of cross-correlations is roughly tested
#    def test_cor_data_symm1(self):
#	for k in Data_dict.keys():
#	    for i in range(3):
#                for j in range(3):
#		    thcor = Data_dict[k].hcor(U_dict[i],U_dict[j])
#		    thcor_symm = Data_dict[k].hcor_symm(U_dict[i],U_dict[j])
#		    self.assertTrue(np.allclose(thcor[1][len(thcor_symm[1])-1:],np.flipud(thcor_symm[1]),atol=1),
#                            "direction = %d, %d" %(i,j))


#TODO - write tests for kinetic energy

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
        self.assertTrue(np.allclose(hs.symm("symm",np.full(14,3.14)),np.full(7,3.14)))

    def test_symm3(self):
        self.assertTrue(np.allclose(hs.symm("asymm",np.full(32,3.1415)),np.zeros(16)))


    def test_utau(self):
        self.assertAlmostEqual(test_128.utau(),0.04285714285714286)


#******************************************
#call of the program
#******************************************

if __name__=='__main__':
    unittest.main()

