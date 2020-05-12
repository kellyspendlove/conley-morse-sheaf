### VertexFlowSheaf.py
### MIT LICENSE 2019 Kelly Spendlove
### Construct assignment on vertices of flow complex given transversality specs

from CubicalSheaf import *

from pychomp._chomp import *
from pychomp.GeometricCubicalComplex import *
from pychomp.TransversalityComplex import *
from pychomp.FlowComplex import *

def TransversalitySheaf(base_complex, ODE_model,phase_specs,num_samples):
  """
  Input:
    complex : Geometric cubical complex as base space (parameter space)
    ODE_model : Right hand side of ODE
    phase_specs : bounds on phase space, number of boxes across
    num_samples : number of samples for codimension 1 walls
  """
  mapping = {}
  #Step 1. Construct the complex for phase space
  bounds,boxes = phase_specs
  geometric_phase_complex = GeometricCubicalComplex(bounds,boxes)
  phase_complex = geometric_phase_complex.complex
  #Step 2. For each cell, compute transversality complex, followed by graded_cell complex
  # 1 . The transversality complex, returns (phase_complex, discrete_flow)
  B = base_complex.complex
  for base_vertex in B(0):
    params = [coord[0] for coord in base_complex.geometry(base_vertex)]
    (C,discrete_flow) = TransversalityComplex(geometric_phase_complex, ODE_model, params, num_samples)
    mapping[base_vertex] = FlowComplex(C,discrete_flow)
    #sheaf_mapping[base_vertex] = FlowGradedComplex(*TransversalityComplex(geometric_phase_complex, ODE_model, params, num_samples))
  return CubicalSheaf(base_complex, mapping) 

