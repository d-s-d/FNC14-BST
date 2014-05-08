#!/usr/bin/python

# Ideally, this script will be used to generate plots in the future

import json
from pprint import pprint

with open('test.log') as json_data:
	data = json.load(json_data)
	pprint(data)
