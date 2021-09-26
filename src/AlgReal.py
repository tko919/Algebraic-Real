from src.Template import *
from src.StrumSeq import *
from src.Factorize import *

def Resultant(f,g):
   n,m=f.deg(),g.deg()
   zero,one=f.id1(),f.id2()
   if g.v==UniPoly:
      zero,one=g.id1(),g.id2()
   L=n+m
   Syl=[[zero]*L for i in range(L)]
   for i in range(m):
      for j in range(n+1):
         Syl[i][i+j]=copy.copy(f[j])
   for i in range(n):
      for j in range(m+1):
         Syl[i+m][i+j]=copy.copy(g[j])
   dp=[[zero]*L for i in range(L)]   
   ret=zero
   for i in range(L):
      dp[i][i]=one
   for Len in range(1,L+1):
      nxt=[[zero]*L for i in range(L)]
      for h in range(L):
         for c in range(h,L):
            if dp[h][c]==zero:
               continue
            for to in range(h+1,L):
               nxt[h][to]+=dp[h][c]*-Syl[c][to]
            if Syl[c][h]==zero:
               continue
            add=dp[h][c]*Syl[c][h]
            if Len==L:
               ret+=add
            for to in range(h+1,L):
               nxt[to][to]+=add
      dp=nxt
   return ret

class AlgReal:
   def __init__(self,f=UniPoly([]),a=-1,b=1,x=0):
      self.f,self.a,self.b,self.x=f,a,b,Fraction(x)
      if not isinstance(f,UniPoly):
         if isinstance(f,AlgReal):
            self.f,self.a,self.b,self.x=f.f,f.a,f.b,f.x
         else:
            self.x=Fraction(f)
            self.f=UniPoly([])
         return
      if not self.f.cs:
         return
      if self.f.value(0)==0:
         self.f=UniPoly([])
         self.x=Fraction(0)
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
         self.x=Fraction(self.b)
      
   def __copy__(self):
      return AlgReal(copy.copy(self.f),self.a,self.b,self.x)
   
   def __str__(self):
      ret='{'
      if self.f.cs:
         ret+=str(self.f)+','
         ret+=str(self.a)+','
         ret+=str(self.b)
      else:
         ret+=str(self.x)
      ret+='}'
      return ret

   def __repr__(self):
      return str(self)

   def __eq__(self,other):
      if not isinstance(other,AlgReal):
         if not self.f.cs:
            return self.x==other
         else:
            return self.f.value(other)==0
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
         if other.a>=self.x:
            return True
         if other.b<self.x:
            return False
         other.f.value(self.x)*other.f.value(other.b)<0
      elif not other.f.cs:
         if self.a>=other.x:
            return False
         if self.b<other.x:
            return True
         self.f.value(other.x)*self.f.value(self.b)>0
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
         for i in range(ret.f.deg()+1):
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
         for i in range(ret.f.deg()+1):
            ret.f[i]*=other.denominator**(ret.f.deg()-i)
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
         ret.f=Resultant(f,g)
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
      if ret.f.cs:
         ret.f=ret.f.squarefree()
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
         ret.f=Resultant(f,g)
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
      if ret.f.cs:
         ret.f=ret.f.squarefree()
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
         seq=StrumSeq(ret.f)
         while True:
            if CountRealroots(ret.f,ret.a,ret.b,seq)==1:
               break
            mid=(ret.a+ret.b)/2
            if self.f.value(ret.a)*self.f.value(mid)<=0:
               ret.b=mid
            else:
               ret.a=mid
         ret.a,ret.b=ret.a**k,ret.b**k
         if ret.a>ret.b:
            ret.a,ret.b=ret.b,ret.a
         for g in Factorize(ret.f):
            if g.value(ret.a)-g.value(ret.b)<=0:
               ret.f=g
               break
         if ret.f.deg()==1:
            ret.x=Fraction(-ret.f[0],ret.f[1])
            ret.f.cs=[]
      return ret

def BoundOfRoots(f):
   if not f.cs:
      return 0
   mx=int(max(np.abs(f.cs[:-1])))
   return Fraction(mx,abs(f.lc()))+1

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

def EnumRootsOfAlgRealCoefficient(f):
   def y_i(f,i):
      if i==-1:
         return f
      else:
         return UniPoly([y_i(f,i-1)])
   n=f.deg()
   F=UniPoly([0])
   for i in range(n+1):
      F=UniPoly([F,y_i(UniPoly([1]).shift(i),i)])
   for i in range(n,-1,-1):
      G=f[i].f
      if not f[i].f.cs:
         G=UniPoly([-f[i].x.numerator,f[i].x.denominator])
      F=Resultant(F,G)
   ret=[]
   g=F[0].squarefree()
   for x in EnumRoots(g):
      if f.value(AlgReal(x.a))*f.value(AlgReal(x.b))<=AlgReal(0):
         ret.append(x)
   return ret