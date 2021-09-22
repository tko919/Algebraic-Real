from src.Template import *
from src.UniPoly import *
from src.GaloisField import *

def CantorZassenhaus(f):
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
      if cur.deg()<d*2:
         if cur.deg():
            ret.append(cur)
         break
      base=base.powmod(GF.p,cur)
      sub=base-UniPoly([0,1],f.v)
      g=gcdP(sub,cur)
      for C in EDF(g,d):
         ret.append(C)
      cur/=g
      d+=1
   return ret

def FactorCoefficientBound(f):
   n=f.deg()
   return 2**n*int(np.sqrt(n+1)+1)*int(np.max(f.cs))

def Factorize(f):
   B=FactorCoefficientBound(f)*abs(f.lc())
   assert(B<2**63)
   f_=UniPoly([],GF)
   while True:
      p=random.randint(B*2+1,(B*2+1)*2)
      for d in range(2,int(np.sqrt(p))+1):
         if p%d==0:
            p=-1
            break
      if p==-1:
         continue
      GF.set_order(p)
      g=UniPoly([],GF)
      for c in f.cs:
         g.cs.append(GF(c))
      if gcdP(g,g.diff()).deg():
         continue
      f_=g
      break
   sub=CantorZassenhaus(f_)
   ret=[]
   S=[]
   for i in range(len(sub)):
      S.append(i)
   k=1
   while True:
      if k*2>len(S):
         ret.append(f)
         break
      for p in itertools.combinations(S,k):
         gsum,hsum=0,0
         g_=UniPoly([1],GF)
         for i in p:
            g_*=sub[i]
         g=UniPoly([])
         for c in g_.cs:
            if c.val>B:
               g.cs.append(c.val-GF.p)
            else:
               g.cs.append(c.val)
            g.cs[-1]*=f.lc()
            gsum+=abs(g.lc())
         h_=UniPoly([1],GF)
         for i in S:
            if not i in p:
               h_*=sub[i]
         h=UniPoly([])
         for c in h_.cs:
            if c.val>B:
               h.cs.append(c.val-GF.p)
            else:
               h.cs.append(c.val)
            h.cs[-1]*=f.lc()
            hsum+=abs(h.lc())
         if gsum*hsum<B:
            for i in p:
               S.remove(i)
            f=h.to_pp()
            ret.append(g.to_pp())
            k-=1
            break
      k+=1
   return ret