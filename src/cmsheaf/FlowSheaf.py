### FlowSheaf.py
### MIT LICENSE 2019 Kelly Spendlove
### Construct sheaf based on flow complexes over vertices

from CubicalSheaf import *

from pychomp._chomp import *
#from pychomp.FlowGradedComplex import *
from pychomp.FlowComplex import *
from pychomp.GradedComplex import *

def SheafFromVertices(vertex_sheaf):
  """
  Input:
    complex : Geometric cubical complex as base space (parameter space)
    ODE_model : Right hand side of ODE
    phase_specs : bounds on phase space, number of boxes across
    num_samples : number of samples for codimension 1 walls
  """
  mapping = {}
  B = vertex_sheaf.base_complex.complex
  D = B.dimension()
  for base_vertex in B(0):
    mapping[base_vertex] = vertex_sheaf.mapping[base_vertex]
  C = mapping[next(B(0))].complex
  vertices = [ cell for cell in C(C.dimension())]
  for d in range(1,D+1):
    for base_cell in B(d):
      boundary_flows = []
      for base_boundary in B.boundary({base_cell}):
        flow_complex = mapping[base_boundary]
        boundary_flows . append (flow_complex.discrete_flow)
      mapping[base_cell] = FlowComplex(C, dict([(v,set.union(*[F[v] for F in boundary_flows])) for v in vertices]))
  return CubicalSheaf(vertex_sheaf.base_complex, mapping) 

#Turn sheaf of flow complexes into sheaf of graded complexes
def GradeFlowSheaf (flow_sheaf):
  mapping = {}
  #B = flow_sheaf.base_complex.complex
  face_poset = flow_sheaf.base_complex.face_poset
  for base_cell in face_poset.vertices():
    #Replace with constructing graded complex induced from flow
    flow_complex = flow_sheaf.mapping[base_cell]
    mapping[base_cell] = GradedComplexObj.induce_from_flow(flow_complex.complex, lambda x: flow_complex.discrete_flow[x])
  for base_edge in face_poset.edges():
    [u,v] = base_edge
    gc_u, gc_v = mapping[u], mapping[v]
    mapping[base_edge] = GradedComplexMap.induce_from_value(gc_v, gc_u)
   #mapping[base_cell] = FlowGradedComplex(*flow_sheaf.mapping(base_cell))
  return CubicalSheaf(flow_sheaf.base_complex, mapping)
