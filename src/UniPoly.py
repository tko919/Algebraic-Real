from src.Template import *

def gcd(x,y):
   return (x if y==0 else gcd(y,x%y))

def SubresultantPRS(f,g):
   ret=[copy.copy(f),copy.copy(g)]
   psi=f.v(-1)
   while True:
      f,g=copy.copy(ret[-2]),copy.copy(ret[-1])
      ff=f.scale(g.lc()**(f.deg()-g.deg()+1))
      add=ff%g
      if not add.cs:
         break
      c=1
      if len(ret)==2:
         if (f.deg()-g.deg()+1)&1:
            c=-1
      else:
         c=-f.lc()*psi**(f.deg()-g.deg())
      add=add.unscale(c)
      ret.append(add)
      if f.v==int:
         psi=(-g.lc())**(f.deg()-g.deg())//(psi**(f.deg()-g.deg()-1))
      else:
         psi=(-g.lc())**(f.deg()-g.deg())/(psi**(f.deg()-g.deg()-1))
   return ret

def gcdP(f,g):
   if f.deg()<g.deg():
      f,g=g,f
   ff=f.scale(g.lc()**(f.deg()-g.deg()+1))
   ff%=g
   if not ff.cs:
      return g.to_pp()
   else:
      ps=SubresultantPRS(g,ff)
      return ps[-1].to_pp()

class UniPoly:
   def __init__(self,cs,Value=int) -> None:
      self.cs=[]
      self.v=Value
      if not isinstance(cs,list):
         cs=[cs]
      for c in cs:
         if isinstance(c,list):
            c=UniPoly(c)
         if isinstance(c,UniPoly):
            self.cs.append(c)
            self.v=UniPoly
         else:
            self.cs.append(Value(c))
      while self.cs and self.lc()==self.id1():
         self.cs.pop(-1)
       
   def __copy__(self):
      return UniPoly(self.cs,self.v)

   def __eq__(self,other):
      if not isinstance(other,UniPoly) or self.deg()!=other.deg():
         return False
      for i in range(self.deg()+1):
         if self.cs[i]!=other.cs[i]:
            return False
      return True

   def lc(self):
      return self.cs[-1]

   def id1(self):
      if self.v==UniPoly:
         return UniPoly([])
      else:
         return self.v(0)
         
   def id2(self):
      if self.v==UniPoly:
         return UniPoly([self.lc().id2()])
      else:
         return self.v(1)
         
   def deg(self):
      return len(self.cs)-1

   def scale(self,a):
      ret=copy.copy(self)
      for i in range(len(self.cs)):
         ret.cs[i]*=a
      return ret

   def unscale(self,a):
      ret=copy.copy(self)
      for i in range(len(self.cs)):
         if ret.v==int:
            ret.cs[i]//=a
         else:
            ret.cs[i]/=a
      return ret

   def to_monic(self):
      if not self.cs:
         return self
      return self.unscale(abs(self.lc()))
   
   def to_pp(self):
      if not self.cs:
         return self
      g=0
      for c in self.cs:
         g=gcd(g,c)
      return self.unscale(abs(g))

   def shift(self,a):
      ret=copy.copy(self)
      ret.cs[0:0]=[ret.id1()]*a
      return ret

   def value(self,x):
      ret=x-x
      for c in reversed(self.cs):
         ret*=x
         ret+=c
      return ret

   def __neg__(self):
      ret=copy.copy(self)
      for i in range(len(self.cs)):
         ret.cs[i]=-ret.cs[i]
      return ret

   def __add__(self,other):
      if not isinstance(other,UniPoly):
         ret=copy.copy(self)
         if not ret.cs:
            ret.cs.append(self.v(0))
         ret.cs[0]+=other
         return ret
      L=max(len(self.cs),len(other.cs))
      ret=[self.id1()]*L
      for i in range(L):
         if i<len(self.cs):
            ret[i]+=self.cs[i]
         if i<len(other.cs):
            ret[i]+=other.cs[i]
      return UniPoly(ret,self.v)

   def __sub__(self,other):
      L=max(len(self.cs),len(other.cs))
      ret=[self.id1()]*L
      for i in range(L):
         if i<len(self.cs):
            ret[i]+=self.cs[i]
         if i<len(other.cs):
            ret[i]-=other.cs[i]
      return UniPoly(ret,self.v)

   def __mul__(self,other):
      if not isinstance(other,UniPoly):
         return self.scale(other)
      if not self.cs or not other.cs:
         return UniPoly([])
      L=len(self.cs)+len(other.cs)-1
      ret=[other.id1()]*L
      for i in range(len(self.cs)):
         for j in range(len(other.cs)):
            if isinstance(self.cs[i],UniPoly):
               ret[i+j]+=self.cs[i]*other.cs[j]
            else:
               ret[i+j]+=other.cs[j]*self.cs[i]
      return UniPoly(ret,self.v)

   def __truediv__(self,other):
      if not self.cs or not other.cs:
         return UniPoly([])
      ret=[self.id1()]*(len(self.cs)-len(other.cs)+1)
      f,g=copy.copy(self),copy.copy(other)
      while len(f.cs)>=len(g.cs):
         p,q=f.lc(),g.lc()
         shift=len(f.cs)-len(g.cs)
         if isinstance(p,int):
            ret[shift]=p//q
         else:
            ret[shift]=p/q
         sub=g.scale(ret[shift]).shift(shift)
         f-=sub
      return UniPoly(ret,self.v)

   def __mod__(self,other):
      h=self/other
      ret=self-h*other
      return ret

   def __pow__(self,k):
      ret=UniPoly([self.id2()])
      b=copy.copy(self)
      while k>0:
         if k&1:
            ret*=b
         b*=b
         k>>=1
      return ret

   def comp(f,g):
      ret=UniPoly([f.id1()])
      for c in reversed(f.cs):
         ret*=g
         ret+=c
      return ret

   def diff(self):
      ret=copy.copy(self)
      if ret.cs:
         ret.cs.pop(0)
      for i in range(len(ret.cs)):
         ret.cs[i]*=(i+1)
      return ret

   def squarefree(self):
      g=self.diff()
      ret=self/gcdP(self,g)
      return ret

def debug(f,End='\n'):
   if isinstance(f,UniPoly):
      print('[',end='')
      for i in range(len(f.cs)):
         debug(f.cs[i],(',' if i!=len(f.cs)-1 else ''))
      print(end=']'+End)
   elif isinstance(f,list):
      print('[',end='')
      for i in range(len(f)):
         debug(f[i],(',' if i!=len(f)-1 else ''))
      print(end=']'+End)
   else:
      print(f,end=End)