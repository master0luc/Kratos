from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of a scalar variable using the BofangConditionTemperatureProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeListOfNodalYoungModulusProcess(Model, settings["Parameters"])

class ImposeListOfNodalYoungModulusProcess(Process):
    
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        settings.AddEmptyValue("is_fixed").SetBool(True)

        # Creating the input table according the input values
        min_young = settings["young_Modulus_min"].GetDouble()
        max_young = settings["young_Modulus_max"].GetDouble()
        number_of_cases = settings["number_of_cases"].GetInt()
        delta = (max_young-min_young)/number_of_cases

        self.table = PiecewiseLinearTable()
        for i in range (number_of_cases):
            self.table.AddRow(float(i+1), float(min_young +(i*delta)))
        self.table.AddRow(float(i+2), float(max_young))

        self.process = DamListTableNodalYoungModulusProcess(model_part, self.table, settings) 
                 
                 
    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()