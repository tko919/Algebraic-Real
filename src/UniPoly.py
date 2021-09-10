from src.Template import *

class UniPoly:
   def __init__(self,cs) -> None:
      self.cs=[]
      for c in cs:
         self.cs.append(Fraction(c))
      while len(self.cs)>0 and self.cs[-1]==0:
         self.cs.pop(-1)
      pass
   def id():
      return UniPoly([0,1])
   def deg(self):
      return len(self.cs)-1
   def scale(self,a):
      ret=UniPoly(self.cs)
      for i in range(len(self.cs)):
         ret.cs[i]*=a
      return ret
   def to_monic(self):
      if len(self.cs)==0:
         return self
      lc = abs(self.cs[-1])
      return UniPoly.scale(self,1/lc)
   def shift(self,a):
      ret=UniPoly(self.cs)
      ret.cs[0:0]=[0]*a
      return ret
   def value(self,x):
      ret=0
      for c in reversed(self.cs):
         ret*=x
         ret+=c
      return ret
   def debug(self):
      print(self.cs)
   def __neg__(self):
      ret=UniPoly(self.cs)
      for i in range(len(self.cs)):
         ret.cs[i]=-ret.cs[i]
      return ret

   def __add__(self,other):
      L=max(len(self.cs),len(other.cs))
      ret=[0]*L
      for i in range(L):
         if i<len(self.cs):
            ret[i]+=self.cs[i]
         if i<len(other.cs):
            ret[i]+=other.cs[i]
      return UniPoly(ret)
   def __sub__(self,other):
      L=max(len(self.cs),len(other.cs))
      ret=[0]*L
      for i in range(L):
         if i<len(self.cs):
            ret[i]+=self.cs[i]
         if i<len(other.cs):
            ret[i]-=other.cs[i]
      return UniPoly(ret)
   def __mul__(self,other):
      L=len(self.cs)+len(other.cs)-1
      ret=[0]*L
      for i in range(len(self.cs)):
         for j in range(len(other.cs)):
            ret[i+j]+=self.cs[i]*other.cs[j]
      return UniPoly(ret)
   def __truediv__(self,other):
      ret=[0]*(len(self.cs)-len(other.cs)+1)
      f,g=UniPoly(self.cs),UniPoly(other.cs)
      while len(f.cs)>=len(g.cs):
         p,q=f.cs[-1],g.cs[-1]
         shift=len(f.cs)-len(g.cs)
         ret[shift]=p/q
         sub=UniPoly(g.cs).scale(p/q).shift(shift)
         f-=sub
      return UniPoly(ret)
   def __mod__(self,other):
      h=self/other
      ret=self-h*other
      return ret
   def comp(f,g):
      ret=UniPoly([0])
      for c in reversed(f.cs):
         ret*=g
         ret=UniPoly.scale(ret,c)
      return ret
   def gcd(f,g):
      if len(g.cs)==0:
         return f
      else:
         h=f%g
         return UniPoly.gcd(g,h.to_monic())
   def diff(self):
      ret=UniPoly(self.cs)
      if len(ret.cs)>0:
         ret.cs.pop(0)
      for i in range(len(ret.cs)):
         ret.cs[i]*=(i+1)
      return ret
   def squarefree(self):
      g=UniPoly.diff(self)
      return self/UniPoly.gcd(self,g)