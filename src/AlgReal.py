from src.Template import *
from src.StrumSeq import *

def Resultant(f:UniPoly,g:UniPoly):
   if not g.cs:
      return f.id1()
   if f.deg()<g.deg():
      ret=Resultant(g,f)
      if (f.deg()*g.deg())&1:
         ret=-ret
      return ret
   k=0
   if isinstance(g.lc(),UniPoly):
      k=f.deg()-g.deg()+1
      f=f.scale(g.lc()**k)
   r=f%g
   if not r.cs:
      if g.deg()==0:
         return g.cs[-1]**f.deg()
      else:
         return f.id1()
   else:
      ret=Resultant(g,r)
      j=f.deg()-r.deg()-g.deg()*k
      if j>=0:
         ret*=g.lc()**j
      else:
         ret/=g.lc()**j
      if (f.deg()*g.deg())&1:
         ret=-ret
      return ret

class AlgReal:
   def __init__(self,f=UniPoly([-1,1]),a=Fraction(1,2),b=Fraction(3,2)) -> None:
      self.f=f
      self.a=a
      self.b=b
      self.x=Inf
      if self.f.value(0)==0:
         self.x=0
         self.f=[]
         return
      while self.f.cs and self.f.cs[0]==0:
         self.f.cs.pop(0)
      while True:
         if 0<self.a or self.b<0:
            break
         mid=(self.a+self.b)
         if self.f.value(self.a)*self.f.value(mid)<0:
            self.b=mid
         else:
            self.a=mid
      if self.f.value(self.b)==0:
         self.x=self.b
         self.f=[]
      
   def __copy__(self):
      return AlgReal(self.f,self.a,self.b)
   
   def __eq__(self,other):
      if self.a>other.b or self.b<other.a:
         return False
      if not self.f and not other.f:
         return self.x==other.x
      elif not self.f:
         return other.f.value(self.x)==0
      elif not other.f:
         return self.f.value(other.x)==0
      return CountRealroots(UniPoly.gcd(self.f,other.f)
      ,max(self.a,other.a),min(self.b,other.b))==1

   def __neg__(self):
      ret=copy.copy(self)
      if not ret.f:
         ret.x=-ret.x
      else:
         for i in range(len(ret.f.cs)):
            if i&1:
               ret.f.cs[i]=-ret.f.cs[i]
         ret.a,ret.b=-ret.b,-ret.a
      return ret

   def __add__(self,other):
      ret=copy.copy(self)
      if not isinstance(other,AlgReal):
         ret.f=UniPoly.comp(ret.f,UniPoly([-other,1]))
         ret.a+=other
         ret.b+=other
      elif not self.f and not other.f:
         ret.x+=other.x
      elif not self.f:
         ret=other+ret.x
      elif not other.f:
         ret+=other.x
      else:
         f=UniPoly.comp(self.f,UniPoly([[0,1],[-1]]))
         g=UniPoly.comp(other.f,UniPoly([[],[1]]))
         ret.f=Resultant(f,g).squarefree()
         n=ret.f.deg()
         D=np.sqrt(float(abs(Resultant(ret.f,ret.f.diff())/ret.f.lc())))
         m=D/(abs(ret.f.lc())**(n-1)*(BoundOfRoots(ret.f)*2)**(n*(n-1)//2-1))
         a,b,c,d=self.a,self.b,other.a,other.b
         while True:
            ret.a,ret.b=a+c,b+d
            if ret.b-ret.a<m:
               break
            mid=(a+b)/2
            if self.f.value(a)*self.f.value(mid)<=0:
               b=mid
            else:
               a=mid
            mid=(c+d)/2
            if other.f.value(c)*other.f.value(mid)<=0:
               d=mid
            else:
               c=mid    
      return ret

   def __sub__(self,other):
      return self+(-other)

   def __mul__(self,other):
      ret=copy.copy(self)
      if not isinstance(other,AlgReal):
         q=1
         for i in range(len(ret.f.cs)):
            ret.f.cs[i]/=q
            q*=other
         ret.a*=other
         ret.b*=other
         if other<0:
            ret.a,ret.b=ret.b,ret.a
      if not self.f and not other.f:
         ret.x*=other.x
      elif not self.f:
         ret=other*ret.x
      elif not other.f:
         ret*=other.x
      else:
         base=[]
         n=self.f.deg()
         for i in range(n+1):
            add=UniPoly([self.f.cs[n-i]])
            base.append(add.shift(n-i))
         f=UniPoly(base)
         g=UniPoly.comp(other.f,UniPoly([[],[1]]))
         ret.f=Resultant(f,g).squarefree()
         n=ret.f.deg()
         D=np.sqrt(float(abs(Resultant(ret.f,ret.f.diff())/ret.f.lc())))
         m=D/(abs(ret.f.lc())**(n-1)*(BoundOfRoots(ret.f)*2)**(n*(n-1)//2-1))
         a,b,c,d=self.a,self.b,other.a,other.b
         while True:
            ret.a,ret.b=min(a*c,a*d,b*c,b*d),max(a*c,a*d,b*c,b*d)
            if ret.b-ret.a<m:
               break
            mid=(a+b)/2
            if self.f.value(a)*self.f.value(mid)<=0:
               b=mid
            else:
               a=mid
            mid=(c+d)/2
            if other.f.value(c)*other.f.value(mid)<=0:
               d=mid
            else:
               c=mid    
      return ret
   
   def inv(self):
      ret=copy.copy(self)
      ret.f.cs.reverse()
      ret.a,ret.b=1/ret.b,1/ret.a
      return ret

   def __truediv__(self,other):
      return self*other.inv()


def EnumRoots(f,a=-Inf,b=Inf,seq=[]):
   if a==-Inf:
      Bound=BoundOfRoots(f)
      a=max(min(a,Bound),-Bound)
      b=max(min(b,Bound),-Bound)
   if not seq:
      seq=StrumSeq(f)
   
   if CountRealroots(f,a,b,seq)==0:
      return []
   elif CountRealroots(f,a,b,seq)==1:
      return [AlgReal(f,a,b)]
   else:
      mid=(a+b)/2
      ret=EnumRoots(f,a,mid,seq)
      for x in EnumRoots(f,mid,b,seq):
         ret.append(x)
      return ret