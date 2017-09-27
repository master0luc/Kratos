# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
# Condition number utility
from test_feast_condition_number_utility import TestConditionNumberUtility as TTestConditionNumberUtility

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    # Exact integration tests
    smallSuite.addTest(TTestConditionNumberUtility('test_cond_matrix_2x2'))
    smallSuite.addTest(TTestConditionNumberUtility('test_cond_matrix_3x3'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            ## SMALL
            TTestConditionNumberUtility,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
