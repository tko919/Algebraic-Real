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
         ret=EDF(h,d)+EDF(g/h,d)
         return ret
      else:
         u=u.powmod((GF.p**d-1)//2,g)
         u-=UniPoly([1],GF)
         h=gcdP(u,g)
         if h.deg() and h!=g:
            ret=EDF(h,d)+EDF(g/h,d)
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
      ret+=EDF(g,d)
      cur/=g
      d+=1
   return ret

def FactorCoefficientBound(f):
   n=f.deg()
   return 2**n*int(np.sqrt(n+1)+1)*int(np.max(f.cs))

def IsPrime(p):
   for d in range(2,round(p**.5)+1):
      if p%d==0:
         return False
   return True

def EnumPrime():
   p=1
   while True:
      p+=2
      if IsPrime(p):
         yield p

def HenselLifting(f,gs,B,p):
   def extgcd(a,b):
      s,xs,ys,t,xt,yt=a,UniPoly([1],GF),UniPoly([0],GF),b,UniPoly([0],GF),UniPoly([1],GF)
      while t.deg()>0:
         q=s/t
         s,xs,ys,t,xt,yt=t,xt,yt,s-q*t,xs-q*xt,ys-q*yt
      return xt/t,yt/t
   GF.set_order(p)
   mid=len(gs)//2
   g,h=UniPoly([1],GF),UniPoly([1],GF)
   for i in range(len(gs)):
      if i<mid:
         g*=gs[i]
      else:
         h*=gs[i]
   g=g.scale(f.lc())
   s,t=extgcd(g,h)
   pl=p
   while True:
      if pl>=B*2+1:
         break
      pl*=pl
      GF.set_order(pl)
      e=UniPoly(f.cs,GF)
      e-=g*h
      q=(s*e)/h
      r=s*e-h*q
      k=g.deg()
      g=(g+t*e+g*q).slice(k+1)
      h+=r
      b=s*g+t*h-1
      c=(s*b)/h
      d=s*b-h*c
      s-=d
      t=(t-t*b-g*c).slice(k+1)
   if len(gs)==1:
      return [h.scale(g.lc())]
   elif len(gs)==2:
      return [g,h]
   else:
      return HenselLifting(g,gs[:mid],B,p)+HenselLifting(h,gs[mid:],B,p)

def Factorize(f):
   B=FactorCoefficientBound(f)*abs(f.lc())
   gen=EnumPrime()
   f_=UniPoly([])
   while True:
      p=next(gen)
      if f.lc()%p==0:
         continue
      GF.set_order(p)
      f_=UniPoly([],GF)
      for c in f.cs:
         f_.cs.append(GF(c))
      f_=f_.to_monic()
      if gcdP(f_,f_.diff()).deg()<=0:
         break
   sub=HenselLifting(f,CantorZassenhaus(f_),B,GF.p)
   k=1
   ret=[]
   T={i for i in range(0,len(sub))}
   while True:
      if k*2>len(T):
         if f.deg()>0:
            ret.append(f)
         break
      for S in itertools.combinations(T,k):
         g_,h_=UniPoly([f.lc()],GF),UniPoly([f.lc()],GF)
         for i in T:
            if i in S:
               g_*=sub[i]
            else:
               h_*=sub[i]
         g,h=UniPoly([]),UniPoly([])
         gsum,hsum=0,0
         for c in g_.cs:
            if c.val<=B:
               g.cs.append(c.val)
            else:
               g.cs.append(c.val-GF.p)
            gsum+=abs(g.lc())
         for c in h_.cs:
            if c.val<=B:
               h.cs.append(c.val)
            else:
               h.cs.append(c.val-GF.p)
            hsum+=abs(h.lc())
         if gsum*hsum<B:
            ret.append(g.to_pp())
            f=h.to_pp()
            T-=set(S)
      k+=1
   return ret