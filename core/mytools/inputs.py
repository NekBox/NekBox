#!/usr/bin/python

import json
from sys import argv

with open(argv[1], 'r') as f:
  params = json.load(f)

values = dict(zip(params.keys(), [pair[0] for pair in params.values()]))
comments = dict(zip([key + "_c" for key in params.keys()], [pair[1] for pair in params.values()]))

everything = dict(values.items() + comments.items())

with open("SIZE_template", "w") as f:
  skeys = sorted(values.keys());
  for key in skeys:
    if not key[0:2] == "__": 
      f.write("  parameter ( {0:<6} = {{{0}!s:4}} ) ! {{{0}_c}} \n".format(key))

with open("SIZE_template", 'r') as f:
  size_template = f.read()
size_contents = size_template.format(**everything)
with open(argv[2], 'w') as f:
  f.write(size_contents)

'''
with open(argv[2], 'w') as f:
  for key in values:
    if not key[0:2] == "__": 
      f.write("  parameter( {0!s:<6} = {1!s:4} ) ! {2!s} \n".format(key, values[key], comments[key+"_c"]))

'''

