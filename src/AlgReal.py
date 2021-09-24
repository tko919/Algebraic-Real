from itertools import permutations
from src.Template import *
from src.StrumSeq import *
from src.Factorize import *

def Resultant(f,g):
   n,m=f.deg(),g.deg()
   L=n+m
   Syl=[[f.id1()]*L for i in range(L)]
   for i in range(m):
      for j in range(n+1):
         Syl[i][i+j]=copy.copy(f[j])
   for i in range(n):
      for j in range(m+1):
         Syl[i+m][i+j]=copy.copy(g[j])

   dp=[[f.id1()]*L for i in range(L)]   
   ret=f.id1()
   for i in range(L):
      dp[i][i]=f.id2()
   for Len in range(1,L+1):
      nxt=[[f.id1()]*L for i in range(L)]
      for h in range(L):
         for c in range(h,L):
            if dp[h][c]==f.id1():
               continue
            for to in range(h+1,L):
               nxt[h][to]+=dp[h][c]*-Syl[c][to]
            if Syl[c][h]==f.id1():
               continue
            add=dp[h][c]*Syl[c][h]
            if Len==L:
               ret+=add
            for to in range(h+1,L):
               nxt[to][to]+=add
      dp=nxt
   return ret

class AlgReal:
   def __init__(self,f=UniPoly([]),a=-1,b=1,x=0) -> None:
      self.f=f
      self.a=a
      self.b=b
      self.x=x
      if self.x!=Inf:
         return
      if self.f.value(0)==0:
         self.f=UniPoly([])
         self.x=0
         return
      while self.f.cs and self.f[0]==0:
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
         self.f=UniPoly([])
         self.x=self.b
   
   def rational(r):
      ret=AlgReal()
      ret.f=UniPoly([])
      ret.x=Fraction(r,1)
      return ret
      
   def __copy__(self):
      return AlgReal(copy.copy(self.f),self.a,self.b,self.x)
   
   def __eq__(self,other):
      if not self.f.cs and not other.f.cs:
         return self.x==other.x
      elif not self.f.cs:
         return other.f.value(self.x)==0
      elif not other.f.cs:
         return self.f.value(other.x)==0
      if self.a>other.b or self.b<other.a:
         return False
      return CountRealroots(gcdP(self.f,other.f)
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
               ret.f[i]=-ret.f[i]
         ret.a,ret.b=-ret.b,-ret.a
      ret.x=-ret.x
      return ret

   def __add__(self,other):
      ret=copy.copy(self)
      if isinstance(other,int):
         ret.f=UniPoly.comp(ret.f,UniPoly([-other,1]))
         ret.a+=other
         ret.b+=other
      elif isinstance(other,Fraction):
         ret.f=UniPoly.comp(ret.f,UniPoly([-other.numerator,other.denominator]))
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
         seq=StrumSeq(ret.f)
         a,b,c,d=self.a,self.b,other.a,other.b
         while True:
            ret.a,ret.b=a+c,b+d
            if CountRealroots(ret.f,ret.a,ret.b,seq)==1:
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
         for g in Factorize(ret.f):
            if g.value(ret.a)*g.value(ret.b)<=0:
               ret.f=g
               break
         if ret.f.deg()==1:
            ret.x=Fraction(-ret.f[0],ret.f[1])
            ret.f.cs=[]
      return ret

   def __sub__(self,other):
      return self+(-other)

   def __mul__(self,other):
      ret=copy.copy(self)
      if isinstance(other,int):
         for i in range(ret.f.deg()+1):
            ret.f[i]*=other**(ret.f.deg()-i)
         ret.a*=other
         ret.b*=other
         if other<0:
            ret.a,ret.b=ret.b,ret.a
      elif isinstance(other,Fraction):
         for i in range(ret.f.deg()+1):
            ret.f[i]*=other.denominator**i*other.numerator**(ret.f.deg()-i)
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
            add=UniPoly([self.f[n-i]])
            base.append(add.shift(n-i))
         f=UniPoly(base)
         g=UniPoly.comp(other.f,UniPoly([[],[1]]))
         ret.f=Resultant(f,g).squarefree()
         seq=StrumSeq(ret.f)
         a,b,c,d=self.a,self.b,other.a,other.b
         while True:
            ret.a,ret.b=min(a*c,a*d,b*c,b*d),max(a*c,a*d,b*c,b*d)
            if CountRealroots(ret.f,ret.a,ret.b,seq)==1:
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
         for g in Factorize(ret.f):
            if g.value(ret.a)*g.value(ret.b)<=0:
               ret.f=g
               break
         if ret.f.deg()==1:
            ret.x=Fraction(-ret.f[0],ret.f[1])
            ret.f.cs=[]
      return ret

   def __truediv__(self,other):
      if isinstance(other,Fraction):
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
         base=base.scale(ret.f.lc()**(base.deg()-ret.f.deg()-1))
         base%=ret.f
         ret=base.value(ret)
      return ret

   def kthRoot(self,k):
      ret=copy.copy(self)
      if not self.f.cs:
         ret.f=UniPoly([1]).shift(k)
         ret.f[0]-=ret.x
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
         for g in Factorize(ret.f):
            if g.value(ret.a)-g.value(ret.b)<=0:
               ret.f=g
               break
         if ret.f.deg()==1:
            ret.x=Fraction(-ret.f[0],ret.f[1])
            ret.f.cs=[]
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
      ret=EnumRoots(f,a,mid,seq)+EnumRoots(f,mid,b,seq)
      return ret