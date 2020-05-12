### CubicalSheaf.py
### MIT LICENSE 2019 Kelly Spendlove

"""
Cubical Sheaf Class
"""

from pychomp._chomp import *


class CubicalSheaf:
  def __init__(self, base_complex, mapping):
    """
    Inputs:
      base_complex : A geometric cubical complex on the base space (parameter space)
      mapping : a dictionary which maps cells to graded complexes
    """
    self.base_complex = base_complex
    self.mapping = mapping

  # @classmethod
  # def induce_from_vertex (cls, source, target): 
  # 	#TODO induce sheaf from data on vertices
