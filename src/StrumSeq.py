from src.Template import *
from src.UniPoly import *

def StrumSeq(f):
   pre=copy.copy(f).squarefree()
   ret=SubresultantPRS(pre,pre.diff())
   for i in range(len(ret)):
      if (i&3)<2:
         ret[i]=-ret[i]
   return ret

def variance(seq,x):
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

def CountRealroots(f,a,b,seq=[]):
   if not seq:
      seq=StrumSeq(f)
   return variance(seq,a)-variance(seq,b)

def approximate(f,a,b,err):
   seq=StrumSeq(f)
   assert(CountRealroots(f,a,b,seq)==1)
   while b-a>=err:
      mid=(a+b)/2
      if f.value(a)*f.value(mid)<=0:
         b=mid
      else:
         a=mid
   return a