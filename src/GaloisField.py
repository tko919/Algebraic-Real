from src.Template import *

class GF:
   p=3
   val=0
   def __init__(self,x):
      if isinstance(x,GF):
         self.val=x.val%GF.p
      else:
         self.val=x%GF.p
      if self.val<0:
         self.val+=GF.p

   def set_order(q):
      GF.p=q
   
   def __eq__(self,other):
      if isinstance(other,GF):
         return self.val==other.val
      else:
         return self.val==other

   def __neg__(self):
      return GF(GF.p-self.val)

   def __add__(self,other):
      if isinstance(other,GF):
         return GF(self.val+other.val)
      else:
         return GF(self.val+other)

   def __sub__(self,other):
      if isinstance(other,GF):
         return GF(self.val-other.val+GF.p)
      else:
         return GF(self.val+other+GF.p)
   
   def __mul__(self,other):
      if isinstance(other,GF):
         return GF(self.val*other.val)
      else:
         return GF(self.val*other)

   def __truediv__(self,other):
      if isinstance(other,GF):
         return GF(self.val*pow(other.val,GF.p-2,GF.p))
      else:
         return GF(self.val*pow(other,GF.p-2,GF.p))

   def __pow__(self,k):
      return GF(pow(self.val,k,GF.p))