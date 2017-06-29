# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homstat_tests.py
# Created by: gemusia
# Creation date: 23-06-2017
# Last modified: 28-06-2017 18:46:48
# Purpose: test of module homstat.py
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import homstat as hs
import unittest

class ChannelTest(unittest.TestCase):

    def test_input0(self):
        with self.assertRaises(ValueError):
            hs.Channel([],[],[],0,0,0)

    def test_input1(self):
        with self.assertRaises(ValueError):
            hs.Channel([[[1]]],[[[2]]],[[[3]]],1,1,1)

    def test_input2(self):
        with self.assertRaises(ValueError):
            hs.Channel([[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],1,2,1)

    def test_input3(self):
        with  self.assertRaises(ValueError):
            hs.Channel([[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],2,2,1)

    def test_inputi4(self):
        with self.assertRaises(ValueError):
            hs.Channel([[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]],2,1,1)

class ChannelTest_functions(unittest.TestCase):
    def test_cheb0(self):
        for x in {32,64,128}:
            self.assertEqual(hs.ChebZeros(x+1), np.cos(y*np.pi/N))
  

if __name__=='__main__':
    unittest.main()



'''
tl1 = np.array([])
tl2 = np.array([[]])
tl3 = np.array([[[1],[2]],[3]])
tl4 = np.array([[[1],[2]],[[3],[4]]])
tl5 = np.array([[[1,1,1],[2,2,2]],[[3,3,3],[4,4,4]]])
tl6 = np.array([[[1,1,1.7],[2,2.7,2]],[[3,3.7,3],[4.7,4,4]],[[5,5.7,5],[6.7,6,6.1]]])

'''

