from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

import math

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeForceAdjointModelProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ImposeForceAdjointModelProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        self.components_process_list = []
        
        self.factor = settings["modulus"].GetDouble();
        self.direction = [settings["direction"][0].GetDouble(),settings["direction"][1].GetDouble(),settings["direction"][2].GetDouble()]
        self.value = [self.direction[0]*self.factor,self.direction[1]*self.factor,self.direction[2]*self.factor]
        
        if abs(self.value[0])>1.0e-15:
            x_params = Parameters("{}")
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddValue("mesh_id",settings["mesh_id"])
            x_params.AddValue("position_x",settings["position_x"])
            x_params.AddValue("position_y",settings["position_y"])
            x_params.AddValue("position_z",settings["position_z"])
            x_params.AddEmptyValue("value").SetDouble(self.value[0])
            x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            x_params.AddValue("table",settings["table"])
            self.components_process_list.append(DamApplyForceBySpatialPositionProcess(model_part, x_params))

        if abs(self.value[1])>1.0e-15:
            y_params = Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddValue("mesh_id",settings["mesh_id"])
            y_params.AddValue("position_x",settings["position_x"])
            y_params.AddValue("position_y",settings["position_y"])
            y_params.AddValue("position_z",settings["position_z"])
            y_params.AddEmptyValue("value").SetDouble(self.value[1])
            y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            y_params.AddValue("table",settings["table"])
            self.components_process_list.append(DamApplyForceBySpatialPositionProcess(model_part, y_params))

        if abs(self.value[2])>1.0e-15:
            z_params = Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddValue("mesh_id",settings["mesh_id"])
            z_params.AddValue("position_x",settings["position_x"])
            z_params.AddValue("position_y",settings["position_y"])
            z_params.AddValue("position_z",settings["position_z"])
            z_params.AddEmptyValue("value").SetDouble(self.value[2])
            z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            z_params.AddValue("table",settings["table"])
            self.components_process_list.append(DamApplyForceBySpatialPositionProcess(model_part, z_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
