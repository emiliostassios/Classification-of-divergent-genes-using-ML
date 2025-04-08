#!/usr/bin/env python3

import pandas as pd
import sys

df = pd.read_csv(sys.argv[1],sep = '\t')

class_col = sys.argv[2]

df['CLASS'] = [class_col] * len(df)

df.to_csv(sys.argv[1], index=False, sep = '\t')
