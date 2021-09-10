from src.Template import *
from src.UniPoly import *

def StrumSeq(f):
   pre=UniPoly(f.cs).squarefree().to_monic()
   cur=pre.diff().to_monic()
   ret=[pre]
   while True:
      if len(cur.cs)==0:
         break
      ret.append(cur)
      nxt=-(pre%cur)
      nxt=nxt.to_monic()
      pre=UniPoly(cur.cs)
      cur=UniPoly(nxt.cs)
   return ret

def variance(f,x):
   seq=StrumSeq(f)
   var=[]
   for F in seq:
      y=F.value(x)
      if y==0:
         continue
      if y<0:
         var.append(-1)
      else:
         var.append(1)
   ret=0
   for i in range(len(var)-1):
      if var[i]*var[i+1]<0:
         ret+=1
   return ret

def CountRealroots(f,a,b):
   return variance(f,a)-variance(f,b)

def approximate(f,a,b,err):
   assert(CountRealroots(f,a,b)==1)
   while b-a>=err:
      mid=(a+b)/2
      if CountRealroots(f,a,mid)==1:
         b=mid
      else:
         a=mid
   return float(a)

class AlgReal:
   def __init__(self,f,a,b) -> None:
      self.f=f
      self.a=a
      self.b=b
      pass
   def __eq__(self,other):
      if self.a>other.b or self.b<other.a:
         return False
      return CountRealroots(UniPoly.gcd(self.f,other.f)
         ,max(self.a,other.a),min(self.b,other.b))==1


def EnumRoots(f,a=-Inf,b=Inf):
   mx=0
   for c in f.cs:
      mx=max(mx,abs(c))
   Bound=mx/abs(f.cs[-1])+1
   a=max(min(a,Bound),-Bound)
   b=max(min(b,Bound),-Bound)
   
   if CountRealroots(f,a,b)==0:
      return []
   elif CountRealroots(f,a,b)==1:
      return [AlgReal(f,a,b)]
   else:
      mid=(a+b)/2
      ret=EnumRoots(f,a,mid)
      for x in EnumRoots(f,mid,b):
         ret.append(x)
      return ret