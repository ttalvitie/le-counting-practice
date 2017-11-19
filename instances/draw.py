#!/usr/bin/python3

import os
import subprocess
import sys

for f in sorted(os.listdir("img/")):
	if not f.endswith(".png"):
		continue
	os.remove("img/" + f)

subprocess.call(["g++", "draw.cpp", "-o", "draw", "-std=c++0x", "-O3"])

posets = [f for f in sorted(os.listdir()) if f.endswith(".txt")]

for f in sorted(posets, key=lambda n: (int(n.rsplit("_", 2)[1]), n)):
	sys.stdout.write(f + ":  ")
	sys.stdout.flush()
	instance = f.rsplit(".", 1)[0]
	dot = instance + ".dot"
	png = "img/" + instance + ".png"
	
	sys.stdout.write("producing .dot... ")
	sys.stdout.flush()
	subprocess.call(["./draw", f, dot])
	
	sys.stdout.write("producing .png... ")
	sys.stdout.flush()
	with open(png, "wb") as fp:
		subprocess.call(["dot", dot, "-Tpng"], stdout=fp)
	
	sys.stdout.write("\n")
	os.remove(dot)

os.remove("draw")
