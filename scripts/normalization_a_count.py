#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import sys


df = pd.read_csv(sys.argv[1], sep = '\t')

df = df[df['alignment_count']!=0]
length = {}
normalized = []
for re in SeqIO.parse(sys.argv[2],'fasta'):
	length[re.id] = len(str(re.seq))

for i in df['query']:
	normalized.append(int((df[df['query'] == i]['alignment_count']))/length[i])
df['alignment_count']=normalized

df['query'] = df['query'].astype(str)

df.to_csv(sys.argv[3],index=False,sep = '\t')

