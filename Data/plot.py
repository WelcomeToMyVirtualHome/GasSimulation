#!/usr/bin/python3
import os
from os import path
import sys
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd

def plot(x_label, y_label, file):
	file_path = path.relpath(file)
	height = 5
	width = 20
	opacity = 1
	data = pd.read_table(file, sep=" ")
	X = data[data.columns[0]].values
	Y = data[data.columns[1]].values
	plt.figure(num=None, figsize=(width, height), dpi=80, facecolor='w', edgecolor='k')
	plt.plot(X,Y,color="b", alpha=opacity, label=file)
	plt.title(file)
	plt.ylabel(y_label)
	plt.xlabel(x_label)
	plt.legend()
	plt.tight_layout()
	png_name = file[:file.find(".dat")]
	plt.savefig(png_name + ".png")

P_label = 'P[1u nm^-1 ps^-2]'
T_label = 'T[K]'
t_label = 't[ps]'

plot("tau[ps]", "sigma^2(E)", "E_var.dat")

plot("a", "V", "lattice_constant.dat")

plot(T_label, P_label, "PT.dat")

plot(t_label, P_label, "P_T=500.dat")
plot(t_label, P_label, "P_T=1000.dat")
plot(t_label, P_label, "P_T=1500.dat")
plot(t_label, P_label, "P_T=2000.dat")

plot(t_label, T_label, "T_T=500.dat")
plot(t_label, T_label, "T_T=1000.dat")
plot(t_label, T_label, "T_T=1500.dat")
plot(t_label, T_label, "T_T=2000.dat")