import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import numpy as np
sns.set_context("talk")

H2_1e12 = pd.read_csv("BialyTest_H2_1.00e+12.csv")
H2_5e11 = pd.read_csv("BialyTest_H2_5.00e+11.csv")
H2_1e11 = pd.read_csv("BialyTest_H2_1.00e+11.csv")
H2_5e10 = pd.read_csv("BialyTest_H2_5.00e+10.csv")
H2_2e10 = pd.read_csv("BialyTest_H2_2.00e+10.csv")
H2_1e10 = pd.read_csv("BialyTest_H2_1.00e+10.csv")

dtmax = [1.e10, 2.e10,5.e10, 1.e11, 5.e11, 1.e12]
H2 = []
H2.append(float(H2_1e10['H2']))
H2.append(float(H2_2e10['H2']))
H2.append(float(H2_5e10['H2']))
H2.append(float(H2_1e11['H2']))
H2.append(float(H2_5e11['H2']))
H2.append(float(H2_1e12['H2']))

fig, ax = plt.subplots(figsize=(8,6))
plt.loglog(dtmax, H2, 'o')
plt.xlim(3.e12, 7.e9)  # decreasing time
plt.xlabel(r"$dt_{max}$")
plt.ylabel(r"$y(H_2)$")
plt.tight_layout()
plt.savefig("H2_vs_dtmax_nosubstep.pdf")

with open('yH2_v_Z_UMIST2012.csv') as bialyData:
    lines = bialyData.readlines()
    bialy_Z = []
    bialy_H2 = []
    for line in lines[6:]:
            p = line.split(',')
            bialy_Z.append(float(p[0]))
            bialy_H2.append(float(p[1]))

plt.figure(figsize=(8,6))
plt.loglog(bialy_Z, bialy_H2, label='Bialy', color='black', linestyle=':', lw=2)
bialy_Z = []
plt.loglog(H2_5e10['Z'], H2_5e10['H2'], '.',label='dt = 5E10 s')
plt.loglog(H2_1e11['Z'], H2_1e11['H2'], 'o',markerfacecolor='none', label='dt = 1E11 s')
plt.loglog(H2_5e11['Z'], H2_5e11['H2'], 'o',markersize=20, markerfacecolor='none', label='dt = 5E11 s')
plt.loglog(H2_1e12['Z'], H2_1e12['H2'], 'o',markersize=30, markerfacecolor='none',label='dt = 1E12 s',zorder=10)
plt.loglog(H2_5e12['Z'], H2_5e12['H2'], 'o',markersize=40, markerfacecolor='none', label='dt = 5E12 s')
plt.legend()
plt.xlabel(r"$Z/Z_{\odot}$")
plt.ylabel(r"$y(H_2)$")
plt.tight_layout()
plt.savefig("H2_accuracyComparison.pdf")

with open('yH_v_Z_UMIST2012.csv') as bialyData:
    lines = bialyData.readlines()
    bialy_Z = []
    bialy_H = []
    for line in lines[6:]:
            p = line.split(',')
            bialy_Z.append(float(p[0]))
            bialy_H.append(float(p[1]))

plt.figure(figsize=(8,6))
plt.loglog(bialy_Z, bialy_H, label='Bialy', color='black', linestyle=':', lw=2)
bialy_Z = []
plt.loglog(H2_5e10['Z'], H2_5e10['H'], '.',label='dt = 5E10')
plt.loglog(H2_1e11['Z'], H2_1e11['H'], 'o',markerfacecolor='none', label='dt = 1E11 s')
plt.loglog(H2_5e11['Z'], H2_5e11['H'], 'o',markersize=20, markerfacecolor='none', label='dt = 5E11 s')
plt.loglog(H2_1e12['Z'], H2_1e12['H'], 'o',markersize=30, markerfacecolor='none',label='dt = 1E12 s',zorder=10)
plt.loglog(H2_5e12['Z'], H2_5e12['H'], 'o',markersize=40, markerfacecolor='none', label='dt = 5E12 s')
plt.legend()
plt.xlabel(r"$Z/Z_{\odot}$")
plt.ylabel(r"$y(H)$")
plt.tight_layout()
plt.savefig("H_accuracyComparison.pdf")
