### Model.py
### MIT LICENSE 2019 Kelly Spendlove

"""
ODE Models
"""

#from pychomp._chomp import *
import json
from scipy import integrate
import numpy as np
import math

import matplotlib.pyplot as plt


class Model:
  def __init__(self, model):
    def cusp ( x, params ):
      a,b = params
      return -x**3 + a*x + b
    def double_cusp ( x, params ):
      a,b = params
      def tent(x):
        if x <= 5:
          return x
        else:
          return 10-x
      return -x**3 +tent(a)*x + b
    def saddle_node ( x, params ):
      a,b = params
      return -x**2 + a*x + b/4
    def toggle ( X, params ):
      n,m = 4,4
      a,b = params
      return [ a / ( 1 + X[1]**n) - X[0], b / (1+X[0]**m) - X[1]]
    def reduced_toggle ( X , params ):
      n = 4
      a, = params
      return [ a / ( 1 + X[1]**n) - X[0], a / (1+X[0]**n) - X[1]]
    def relaxation_osc ( X, params):
      n = 2
      c = 5
      a,b=params
      return [a + b * X[1]**n / (1+X[1]**n) - X[0], 
                  c*(a + b * X[1]**n / (1 + X[0]**n+X[1]**n) - X[1] )]
    def repressilator ( X, params ):
      n, b = params
      return [ b / (1+X[2]**n) - X[0], b / (1+X[0]**n) - X[1],
                  b / (1+X[1]**n) - X[2]]
    def replicator1D ( x, params ):
      a,b = params
      return x*(1-x)*(x*(a-b)+b)
    def replicator1DTemp ( x, params ):
      a,b = 5,-5
      T = params
      return x*(1-x)*(x*(a-b)+b)+T*(1-x)
    def replicator2D ( X, params ):
      """
      Replicator, bimatrix game
      """
      a,b = params
      c,d = -2,1
      #A,B = np.array([[0,a],[b,0]]),np.array([[0,c],[d,0]])
      return [X[0]*(1-X[0])*(b+(a-b)*X[1]),
              X[1]*(1-X[1])*(d+(c-d)*X[0])]
    def qtLearning( X, params ):
      #A,B = np.array([[10,0],[0,5]]),np.array([[2,0],[0,4]])
      a,b = -10, -2
      c,d = -2,-15
      Tx,Ty = params
      #print(tempX,tempY)
      return [X[0]*(1-X[0])*(a-(a+b)*X[1]+Tx*(1-X[0])*math.log(abs((1-X[0])/X[0])) ),
              X[1]*(1-X[1])*(c-(c+d)*X[0]+Ty*(1-X[1])*math.log(abs((1-X[1])/X[1])) )]
    def QTLearningFull ( X, params):
      A,B = np.array([[10,0],[0,5]]),np.array([[2,0],[0,4]])
      Tx,Ty = params
      x0,x1,y0,y1 = X
      x,y = np.array([x0,x1]),np.array([y0,y1])
      Bx, Ay = np.matmul(B,x), np.matmul(A,y)
      # return [x0*(Ay[0]-np.dot(x,Ay)),
      #       x1*(Ay[1]-np.dot(x,Ay)),
      #       y0*(Bx[0]-np.dot(y,Bx)),
      #       y1*(Bx[1]-np.dot(y,Bx)) ]
      #print(x0,x1,y0,y1)
      return [x0*(Ay[0]-np.dot(x,Ay)+Tx*x1*math.log(abs(x1/x0))),
              x1*(Ay[1]-np.dot(x,Ay)+Tx*x0*math.log(abs(x0/x1))),
              y0*(Bx[0]-np.dot(y,Bx)+Ty*y1*math.log(abs(y1/y0))),
              y1*(Bx[1]-np.dot(y,Bx)+Ty*y0*math.log(abs(y0/y1))) ]
    if model.lower() == 'cusp':
      self.ODE_Model = lambda x, params: cusp( x, params )
      config_file = 'configs/cusp_config.json'
    elif model.lower() == 'double_cusp':
      self.ODE_Model = lambda x, params: double_cusp( x, params )
      config_file = 'configs/double_cusp_config.json'
    elif model.lower() =='sn':
      self.ODE_Model = lambda x, params: saddle_node( x, params )
      config_file = 'configs/SN_config.json'
    elif model.lower() == 'toggle':
      self.ODE_Model = lambda x, params: toggle ( x, params )
      config_file = 'configs/toggle_config.json'
    elif model.lower() == 'reduced_toggle':
      self.ODE_Model = lambda x, params : reduced_toggle ( x, params)
      config_file = 'configs/reduced_toggle_config.json'
    elif model.lower() == 'relaxation_osc':
      self.ODE_Model = lambda x, params : relaxation_osc ( x, params)
      config_file = 'configs/relaxation_osc_config.json'
    elif model.lower() == 'repressilator':
      self.ODE_Model = lambda x, params : repressilator ( x, params)
      config_file = 'configs/repressilator_config.json'
    elif model.lower() == 'qtlearning':
      self.ODE_Model = lambda x, params : qtLearning (x, params )
      config_file = 'configs/qtlearning.json'
    elif model.lower() == 'qtlearningfull':
      self.ODE_Model = lambda x, params : QTLearningFull (x, params )
      config_file = 'configs/qtlearningFull.json'
    elif model.lower() == 'replicator2d':
      self.ODE_Model = lambda x, params : replicator2D (x, params )
      config_file = 'configs/replicator2D.json'
    elif model.lower() == 'replicator1d':
      self.ODE_Model = lambda x, params : replicator1D (x, params )
      config_file = 'configs/replicator1D.json'
    elif model.lower() == 'replicator1dtemp':
      self.ODE_Model = lambda x, params : replicator1DTemp (x, params )
      config_file = 'configs/replicator1DTemp.json'
    #Load parameters
    with open(config_file) as json_config_file:
      data = json.load(json_config_file)
    self.base_bounds = data['param']['bounds']
    self.base_boxes = data['param']['boxes']
    self.phase_bounds = data['phase']['bounds']
    self.phase_boxes = data['phase']['boxes']
    self.num_samples = data['sampling']['num_samples']

  def simulate (self, init_cond, params):
    # time points
    t = np.linspace(0,20)
    # solve ODE
    y = integrate.odeint(lambda x,t : self.ODE_Model(x,params), init_cond, t )

    # plot results
    plt.plot(t,y)
    plt.xlabel('time')
    plt.ylabel('x(t)')
    plt.show()

 
