from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

class TestConditionNumberUtility(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def test_cond_matrix_2x2(self):

        matrix = KratosMultiphysics.Matrix(2, 2)
        matrix[0, 0] = 0.63672
        matrix[0, 1] = 0.60744  
        matrix[1, 0] = 0.26959  
        matrix[1, 1] = 0.20759  

        cond_util = ExternalSolversApplication.FEASTConditionNumberUtility()
        cond = cond_util.GetConditionNumber(matrix)

        self.assertAlmostEqual(cond, 28.150, 3)
        
    def test_cond_matrix_3x3(self):

        matrix = KratosMultiphysics.Matrix(3, 3)
        matrix[0, 0] = 0.63672
        matrix[0, 1] = 0.60744  
        matrix[0, 2] = 0.48016
        matrix[1, 0] = 0.26959  
        matrix[1, 1] = 0.20759  
        matrix[1, 2] = 0.99380
        matrix[2, 0] = 0.28495  
        matrix[2, 1] = 0.49421  
        matrix[2, 2] = 0.28403

        cond_util = ExternalSolversApplication.FEASTConditionNumberUtility()
        cond = cond_util.GetConditionNumber(matrix)

        self.assertAlmostEqual(cond, 10.967, 3)

if __name__ == '__main__':
    KratosUnittest.main()
                
                
                
                

 
