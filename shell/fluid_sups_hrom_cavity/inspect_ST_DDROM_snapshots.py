#!/usr/bin/env python3
from pathlib import Path
import argparse
p=argparse.ArgumentParser()
p.add_argument('directory', type=Path)
p.add_argument('--ranks', type=int, default=4)
a=p.parse_args()
for rank in range(a.ranks):
    count=0
    while (a.directory/f"fluid_st_window_{count:04d}_rank_{rank:06d}.bin").is_file(): count+=1
    print(f"rank {rank}: consecutive windows={count}")
