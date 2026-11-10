import random, sys, math, pickle
from collections import Counter
from multiprocessing import Pool
r0,s0=int(sys.argv[1]),int(sys.argv[2]); N=int(sys.argv[3]); bits=int(sys.argv[4])
d=2**s0-3**r0
def land(seed):
    rnd=random.Random(seed); out=Counter()
    lim=d*(1<<14)
    for _ in range(N//64):
        x=rnd.getrandbits(bits)|1
        while math.gcd(x,d)!=1: x=rnd.getrandbits(bits)|1
        # descend
        while x>lim:
            x=(3*x+d)>>1 if x&1 else x>>1
        # Brent cycle detection
        power=lam=1; tort=x; x=(3*x+d)>>1 if x&1 else x>>1
        steps=0
        while tort!=x:
            if power==lam: tort=x; power*=2; lam=0
            x=(3*x+d)>>1 if x&1 else x>>1; lam+=1; steps+=1
            if steps>5_000_000: tort=None; break
        if tort is None: out[-1]+=1; continue
        m=x;y=(3*x+d)>>1 if x&1 else x>>1
        while y!=x:
            if y<m: m=y
            y=(3*y+d)>>1 if y&1 else y>>1
        out[m]+=1
    return out
if __name__=="__main__":
    with Pool(64 if False else None) as p: res=p.map(land,range(64))
    tot=Counter()
    for c in res: tot.update(c)
    # cycle stats
    rows=[]
    for m,cnt in tot.items():
        if m<0: rows.append((m,cnt,0,0,0,0)); continue
        x=m;odd=[];evn=0
        while True:
            if x&1: odd.append(x); x=(3*x+d)>>1
            else: evn+=1; x>>=1
            if x==m: break
        h1=sum(d/v for v in odd); h2=sum((d/v)**2 for v in odd)
        rows.append((m,cnt,len(odd),len(odd)+evn,h1,h2))
    pickle.dump((d,rows),open(f'basin_{r0}_{s0}.pkl','wb'))
    print('d',d,'samples',sum(tot.values()),'cycles hit',len(tot),'nonterminating',tot.get(-1,0))
    rows.sort(key=lambda z:-z[1])
    for z in rows[:15]: print(z[1],'r',z[2],'len',z[3],'min/d',round(z[0]/d,3),'h1',round(z[4],4),'h2',round(z[5],5))
