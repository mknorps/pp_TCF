# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: homfigs_tests.py
# Created by: gemusia
# Creation date: 02-07-2017
# Last modified: 02-07-2017 16:08:37
# Purpose: 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
import homfigs as hf
import unittest

ff = hf.Homfig(title='very important title',xlabel="xxx",ylabel="yyy")
ff.hdraw()

class HomfigTest(unittest.TestCase):

    def test_input0(self):
        with self.assertRaises(ValueError):
            hs.Channel([],[],[],0,0,0)
