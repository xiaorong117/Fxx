#!/usr/bin/env python3
"""Regression: immutable AMGX steps may predate step_wall_seconds."""
import csv
rows=[{'assembly_seconds':'1','solve_seconds':'2'}]
fields=['step_wall_seconds','assembly_seconds','solve_seconds']
sums={key:sum(float(row.get(key,0.0)) for row in rows) for key in fields}
assert sums=={'step_wall_seconds':0.0,'assembly_seconds':1.0,'solve_seconds':2.0}
print('PASS optional legacy AMGX step_wall_seconds schema')
