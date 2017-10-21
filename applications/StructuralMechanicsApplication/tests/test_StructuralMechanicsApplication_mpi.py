# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
## SMALL 

## NIGTHLY TESTS
# Shell test

from NightlyTestsMPI import ShellT3ThickLinearStaticTests     as TShellT3ThickLinearStaticTests
from NightlyTestsMPI import ShellT3ThickNonLinearStaticTests  as TShellT3ThickNonLinearStaticTests
from NightlyTestsMPI import ShellT3ThickLinearDynamicTests    as TShellT3ThickLinearDynamicTests
from NightlyTestsMPI import ShellT3ThickNonLinearDynamicTests as TShellT3ThickNonLinearDynamicTests

from NightlyTestsMPI import ShellQ4ThinLinearStaticTests      as TShellQ4ThinLinearStaticTests
from NightlyTestsMPI import ShellQ4ThinNonLinearStaticTests   as TShellQ4ThinNonLinearStaticTests
from NightlyTestsMPI import ShellQ4ThinLinearDynamicTests     as TShellQ4ThinLinearDynamicTests
from NightlyTestsMPI import ShellQ4ThinNonLinearDynamicTests  as TShellQ4ThinNonLinearDynamicTests

## VALIDATION TESTS

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

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    # Shell tests
    nightSuite.addTest(TShellT3ThickLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellT3ThickNonLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellT3ThickLinearDynamicTests('test_execution'))
    nightSuite.addTest(TShellT3ThickNonLinearDynamicTests('test_execution'))

    nightSuite.addTest(TShellQ4ThinLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinNonLinearStaticTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinLinearDynamicTests('test_execution'))
    nightSuite.addTest(TShellQ4ThinNonLinearDynamicTests('test_execution'))
    # CL tests
    ##nightSuite.addTest(TIsotropicDamageSimoJuPSTest('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    # SPRISM tests
    ####validationSuite.addTest(TSprismPanTests('test_execution'))
    
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            # TShellT3ThickLinearStaticTests,
            # TShellT3ThickNonLinearStaticTests,
            TShellT3ThickLinearDynamicTests
            # TShellT3ThickNonLinearDynamicTests,
            # TShellQ4ThinLinearStaticTests,
            # TShellQ4ThinNonLinearStaticTests,
            # TShellQ4ThinLinearDynamicTests,
            # TShellQ4ThinNonLinearDynamicTests
        ])
    )
    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
