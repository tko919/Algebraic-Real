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
      if isinstance(ret,UniPoly):
         ret=ret.to_monic()
      return ret

class AlgReal:
   def __init__(self,f=UniPoly([-1,1]),a=Fraction(1,2),b=Fraction(3,2),x=Inf) -> None:
      self.f=f
      self.a=a
      self.b=b
      self.x=x
      if self.x!=Inf:
         return
      if self.f.value(0)==0:
         self.x=0
         self.f=UniPoly([])
         return
      while self.f.cs and self.f.cs[0]==0:
         self.f.cs.pop(0)
      while True:
         if 0<self.a or self.b<0:
            break
         mid=(self.a+self.b)/2
         if self.f.value(self.a)*self.f.value(mid)<0:
            self.b=mid
         else:
            self.a=mid
      if self.f.value(self.b)==0:
         self.x=self.b
         self.f=UniPoly([])
   
   def rational(r):
      ret=AlgReal()
      ret.f=UniPoly([])
      ret.x=Fraction(r,1)
      return ret
      
   def __copy__(self):
      return AlgReal(self.f,self.a,self.b,self.x)
   
   def __eq__(self,other):
      if not self.f.cs and not other.f.cs:
         return self.x==other.x
      elif not self.f.cs:
         return other.f.value(self.x)==0
      elif not other.f.cs:
         return self.f.value(other.x)==0
      if self.a>other.b or self.b<other.a:
         return False
      return CountRealroots(UniPoly.gcd(self.f,other.f)
      ,max(self.a,other.a),min(self.b,other.b))==1
   
   def __lt__(self,other):
      if not self.f.cs and not other.f.cs:
         return self.x<other.x
      elif not self.f.cs:
         other.f.value(self.x)*other.f.value(other.b)<0
      elif not other.f.cs:
         self.f.value(other.x)*self.f.value(self.b)<0
      else:
         if self==other:
            return False
         if self.b<=other.a:
            return True
         elif other.b<=self.a:
            return False
         a,b=max(self.a,other.a),min(self.b,other.b)
         while True:
            mid=(a+b)/2
            if self.f.value(mid)*self.f.value(self.b)<0:
               if other.f.value(mid)*other.f.value(other.b)<0:
                  a=mid
               else:
                  return False
            else:
               if other.f.value(mid)*other.f.value(other.b)<0:
                  return True
               else:
                  b=mid

   def __gt__(self,other):
      return other<self
   
   def __le__(self,other):
      return self<other or self==other
   
   def __ge__(self,other):
      return other<self or self==other

   def __neg__(self):
      ret=copy.copy(self)
      if ret.f:
         for i in range(len(ret.f.cs)):
            if i&1:
               ret.f.cs[i]=-ret.f.cs[i]
         ret.a,ret.b=-ret.b,-ret.a
      ret.x=-ret.x
      return ret

   def __add__(self,other):
      ret=copy.copy(self)
      if not isinstance(other,AlgReal):
         ret.f=UniPoly.comp(ret.f,UniPoly([-other,1]))
         ret.a+=other
         ret.b+=other
      elif not self.f.cs and not other.f.cs:
         ret.x+=other.x
      elif not self.f.cs:
         ret=other+ret.x
      elif not other.f.cs:
         ret+=other.x
      else:
         f=UniPoly.comp(self.f,UniPoly([[0,1],[-1]]))
         g=UniPoly.comp(other.f,UniPoly([[],[1]]))
         ret.f=Resultant(f,g).squarefree()
         n=ret.f.deg()
         D=abs(Resultant(ret.f,ret.f.diff()))/ret.f.lc()
         m=D/(abs(ret.f.lc())**(n-1)*(BoundOfRoots(ret.f)*2)**(n*(n-1)//2-1))**2
         a,b,c,d=self.a,self.b,other.a,other.b
         while True:
            ret.a,ret.b=a+c,b+d
            if (ret.b-ret.a)**2<m:
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
         q=Fraction(1,1)
         for i in range(len(ret.f.cs)):
            ret.f.cs[i]/=q
            q*=other
         ret.a*=other
         ret.b*=other
         if other<0:
            ret.a,ret.b=ret.b,ret.a
      elif not self.f.cs and not other.f.cs:
         ret.x*=other.x
      elif not self.f.cs:
         ret=other*ret.x
      elif not other.f.cs:
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
         D=abs(Resultant(ret.f,ret.f.diff()))/ret.f.lc()
         m=D/(abs(ret.f.lc())**(n-1)*(BoundOfRoots(ret.f)*2)**(n*(n-1)//2-1))**2
         a,b,c,d=self.a,self.b,other.a,other.b
         while True:
            ret.a,ret.b=min(a*c,a*d,b*c,b*d),max(a*c,a*d,b*c,b*d)
            if (ret.b-ret.a)**2<m:
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

   def __truediv__(self,other):
      if not isinstance(other,AlgReal):
         return self*(1/other)
      elif not other.f.cs:
         base=copy.copy(other)
         base.x=1/base.x
         return self*base
      else:
         base=copy.copy(other)
         base.f.cs.reverse()
         base.a,base.b=1/base.b,1/base.a
         return self*base

   def __pow__(self,k):
      ret=copy.copy(self)
      if not self.f:
         ret.x=ret.x**k
      else:
         base=UniPoly([1]).shift(k)
         base%=ret.f
         ret=base.value(ret)
      return ret

   def kthRoot(self,k):
      ret=copy.copy(self)
      if not self.f.cs:
         ret.f=UniPoly([1]).shift(k)
         ret.f.cs[0]-=ret.x
         if ret.x<0:
            ret.a=min(ret.x,-1)
            ret.b=0
         else:
            ret.a=0
            ret.b=max(ret.x,1)
      else:
         ret.f=UniPoly.comp(ret.f,UniPoly([1]).shift(k))
         for x in EnumRoots(ret.f):
            xk=x**k
            if xk>ret.a and xk<=ret.b:
               ret.a=x.a
               ret.b=x.b
               break    
      return ret


def EnumRoots(f,a=-Inf,b=Inf,seq=[]):
   if a==-Inf:
      Bound=BoundOfRoots(f)
      a=max(min(a,Bound),-Bound)
      b=max(min(b,Bound),-Bound)
   if not seq:
      seq=StrumSeq(f)

   cnt=CountRealroots(f,a,b,seq)
   if cnt==0:
      return []
   elif cnt==1:
      return [AlgReal(f,a,b)]
   else:
      mid=(a+b)/2
      ret=EnumRoots(f,a,mid,seq)
      for x in EnumRoots(f,mid,b,seq):
         ret.append(x)
      return ret