from src.Template import *
from src.UniPoly import *
from src.GaloisField import *

def CantorZassenhaus(f):
   f=f.squarefree()
   #EDF
   def EDF(g,d):
      if g.deg()<d:
         return []
      if g.deg()==d:
         return [g]
      base=[]
      while not base:
         for _ in range(g.deg()+1):
            base.append(GF(random.randint(0,GF.p)))
      u=UniPoly(base,GF)
      h=gcdP(u,g)
      if h.deg():
         ret=EDF(h,d)
         for C in EDF(g/h,d):
            ret.append(C)
         return ret
      else:
         u=u.powmod((GF.p**d-1)//2,g)
         if not u.cs:
            u.cs.append(GF(0))
         u-=UniPoly([1],GF)
         h=gcdP(u,g)
         if h.deg() and h!=g:
            ret=EDF(h,d)
            for C in EDF(g/h,d):
               ret.append(C)
            return ret
         else:
            return EDF(g,d)

   # DDF
   base=UniPoly([0,1],f.v)
   cur=copy.copy(f)
   ret=[]
   d=1
   while True:
      if cur.deg()==0:
         break
      base=base.powmod(GF.p,cur)
      sub=base-UniPoly([0,1],f.v)
      g=gcdP(sub,cur)
      for C in EDF(g,d):
         ret.append(C)
      cur/=g
      d+=1
   return ret