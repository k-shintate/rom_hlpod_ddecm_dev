#!/usr/bin/env python3
import argparse, math
from pathlib import Path
import numpy as np

def read_nodes(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
        return np.array([[float(v) for v in f.readline().split()[:3]] for _ in range(n)])

def read_prisms(path):
    with open(path) as f:
        ne = int(f.readline().split()[0])
        conn = []
        for e in range(ne):
            row = f.readline().split()
            if not row:
                raise RuntimeError(f"truncated {path} at element {e}")
            if len(row) >= 8 and int(row[1]) == 6:
                conn.append([int(v) for v in row[2:8]])
            else:
                conn.append([int(v) for v in row[:6]])
        return np.asarray(conn, dtype=int)

def dshape(r,s,z):
    c=0.5
    dr=np.array([-c*(1-z),c*(1-z),0,-c*(1+z),c*(1+z),0])
    ds=np.array([-c*(1-z),0,c*(1-z),-c*(1+z),0,c*(1+z)])
    dz=np.array([-c*(1-r-s),-c*r,-c*s,c*(1-r-s),c*r,c*s])
    return np.column_stack((dr,ds,dz))

def qpoints():
    a=1/math.sqrt(3)
    pts=[]
    for u in (-a,a):
        for v in (-a,a):
            for w in (-a,a):
                pts.append((0.5*(1+u),0.25*(1-u)*(1+v),w))
    pts.append((1/3,1/3,0.0))
    return pts

def classify(ds,tol):
    mn,mx=float(ds.min()),float(ds.max())
    if mn>tol: return "positive"
    if mx<-tol: return "reversed"
    if mn<-tol and mx>tol: return "sign-changing"
    return "near-zero"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("part_dir",type=Path)
    ap.add_argument("--nranks",type=int,required=True)
    ap.add_argument("--node-base",default="node.dat")
    ap.add_argument("--prism-base",default="graph_elem_prism.dat")
    ap.add_argument("--tol",type=float,default=1e-12)
    ap.add_argument("--show",type=int,default=20)
    args=ap.parse_args()
    for rank in range(args.nranks):
        npth=args.part_dir/f"{args.node_base}.{rank}"
        epth=args.part_dir/f"{args.prism_base}.{rank}"
        if not npth.exists() or not epth.exists():
            print(f"[rank {rank}] missing file"); continue
        xyz=read_nodes(npth); conn=read_prisms(epth)
        if len(conn)==0:
            print(f"[rank {rank}] EMPTY"); continue
        counts={"positive":0,"reversed":0,"sign-changing":0,"near-zero":0}
        bad=[]
        for e,c in enumerate(conn):
            if np.any(c<0) or np.any(c>=len(xyz)):
                bad.append((e,"bad-node-id",float("nan"),float("nan"),c.tolist())); continue
            X=xyz[c]
            ds=np.array([np.linalg.det(X.T@dshape(*q)) for q in qpoints()])
            cls=classify(ds,args.tol); counts[cls]+=1
            if cls!="positive": bad.append((e,cls,float(ds.min()),float(ds.max()),c.tolist()))
        print(f"[rank {rank}] prisms={len(conn)} positive={counts['positive']} reversed={counts['reversed']} sign-changing={counts['sign-changing']} near-zero={counts['near-zero']}")
        for e,cls,mn,mx,c in sorted(bad,key=lambda x:(math.inf if math.isnan(x[2]) else x[2]))[:args.show]:
            print(f"  e={e:6d} {cls:13s} detJ=[{mn:.6e},{mx:.6e}] conn={c}")
        if counts["reversed"]:
            print("  NOTE: reversed -> orientation permutation may fix.")
        if counts["sign-changing"]:
            print("  NOTE: sign-changing -> geometrically folded; reordering alone cannot fix.")
if __name__=="__main__":
    main()

