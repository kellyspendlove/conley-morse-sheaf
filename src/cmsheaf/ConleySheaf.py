### ConleySheaf.py
### MIT LICENSE 2019 Kelly Spendlove

from pychomp._chomp import *

from CubicalSheaf import *
from pychomp.GradedComplex import * 
from pychomp.InducedSubgraph import *
from pychomp.TransitiveClosure import *
from pychomp.Poset import *

# from pychomp.GeometricCubicalComplex import *
# from pychomp.TransversalityComplex import *

def ConleySheaf(cubical_sheaf):
  """
  Input:
    Cubical Sheaf
  Output:
    Conley Sheaf
  """
  conley_mapping = {}
  face_poset = cubical_sheaf.base_complex.face_poset
  for base_cell in face_poset.vertices():
    graded_complex = cubical_sheaf.mapping[base_cell]
    conley_complex = ConnectionMatrix(graded_complex.graded_complex)
    conley_mapping[base_cell] = GradedComplexObj( conley_complex, graded_complex.poset )
  for base_edge in face_poset.edges():
    [u,v] = base_edge
    gc_u, gc_v = cubical_sheaf.mapping[u], cubical_sheaf.mapping[v]
    tower_A, tower_B = ConnectionMatrixTower(gc_u.graded_complex), ConnectionMatrixTower(gc_v.graded_complex)
    conley_mapping[base_edge] = GradedComplexMap.induce_from_tower(tower_B, tower_A, cubical_sheaf.mapping[base_edge].poset_map) 

  return CubicalSheaf(cubical_sheaf.base_complex, conley_mapping) 

def RecurrentSheaf (cubical_sheaf):
  recurrent_mapping = {}
  face_poset = cubical_sheaf.base_complex.face_poset
  for base_cell in face_poset.vertices():
    gc = cubical_sheaf.mapping[base_cell]
    recurrent_poset = Poset(InducedSubgraph(TransitiveClosure(gc.poset.get_children()), lambda v : v in gc.graded_complex.count()))
    recurrent_mapping[base_cell] = GradedComplexObj ( gc.graded_complex, recurrent_poset)
  for base_edge in face_poset.edges():
    recurrent_mapping[base_edge] = cubical_sheaf.mapping[base_edge]

  return CubicalSheaf(cubical_sheaf.base_complex, recurrent_mapping)

    #


    #    rc_u = InducedPoset(gc_u.poset, lambda v : v in gc_u.graded_complex.count() )
