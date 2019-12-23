#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
To run the pgm:
./plot_ICE_spot.py All_data-xxxx.txt
"""


import os
import sys
import matplotlib
matplotlib.rcdefaults()

matplotlib.use('TkAgg')

#['pdf', 'pgf', 'Qt4Agg', 'GTK', 'GTKAgg', 'ps', 'agg', 'cairo', 'MacOSX', 'GTKCairo', 'WXAgg', 'template', 'TkAgg', 'GTK3Cairo', 'GTK3Agg', 'svg', 'WebAgg', 'CocoaAgg', 'emf', 'gdk', 'WX']
import matplotlib.pyplot as plt
import matplotlib.text as mt
import numpy as np
import copy
import re
from matplotlib.widgets import CheckButtons, Button
import pandas as pd
pd.set_option("display.width", 200)

## in case where I start pylab, sum become the sum function from numpy and not the built-in one which is expected.
try:
	del sum
except:
	pass

############### ---------- Functions ---------- ###############
def uniq(seq):
	# Not order preserving
	return {}.fromkeys(seq).keys()


def open_file(file, split_line = "\n", tab = True, sep = "\t"):
	"""Opens a file and returns the content in a list of lines
	If tab=True, we consider the file being a table, so lines are split"""

	f_in = open(file, 'r')

	data = f_in.read().split(split_line)[:-1] # must exist a last blank line
	if tab :
		isfloat = re.compile("^[\+\-]*[\d]+\.?[\d]*$")
		data = [[float(i) if isfloat.match(i.strip()) else i for i in f.split(sep)] for f in data]
	f_in.close()
	return data

def ok_button_callback(event):
	selection=[j for j,i in enumerate(check.lines) if i[0].get_visible()]
	if len(selection)!=0:
		f_out = open("Selection_" + id_locus + ".txt", "a")
		f_out.write(All_data[0]+"\n\n")
#		f_out.write(All_data[1]+"\n\n")
		for s in selection:
			f_out.write(All_data[s+1])
			f_out.write("\n\n")
		f_out.close()
		ok_button.label.set_text("Exported")
	else:
		ok_button.label.set_text("Ok")

def erase_last_print(size):
	CURSOR_UP_ONE = '\x1b[1A'
	ERASE_LINE = '\x1b[2K'
	for i in range(size):
		sys.stdout.write(CURSOR_UP_ONE + ERASE_LINE)

############### ---------- Definitions ---------- ###############


file_in = sys.argv[-1]
id_locus = file_in.split(".")[0].split("-")[-1]
All_data = open_file(file_in, split_line = "\n\n" , tab = 0 )

############### ---------- Read Data ---------- ###############


digit = re.compile("^[\+\-]*[\d]+$")

alldat = All_data[1:]
sizedata = len(All_data[1:])
GC = np.ndarray(sizedata,dtype=np.object)
NAME = np.ndarray(sizedata,dtype=np.object)
STRAIN = np.ndarray(sizedata,dtype=np.object)
dic_data = np.ndarray(sizedata,dtype=np.object)

for j in xrange(sizedata):
	NAME[j] = alldat[j].split("\n")[0]
	STRAIN[j]= {s.split("&")[0]:int(s.split("&")[1]) for s in alldat[j].split("\n")[1].split()}
	dic_data[j] = {data.split("\t")[0] : [int(d) if digit.match(d) else [int(rna) if digit.match(rna) else rna for rna in d.split("&")] if d.find("&") != -1 or d[0:4]==STRAIN[j].keys()[0][0:4] else d.strip() for d in data.split("\t")][1:] for data in alldat[j].split("\n")[2:-1]}
	GC[j] = np.array([float(gc) for gc in alldat[j].split("\n")[-1].split()])

ALL_STRAINS = np.sort(uniq(sum([i.keys() for i in STRAIN],[])))



txt = []
clicks = []
patch_clicked = []
LAST = 0
def onpick1(event):
	global txt
	global count
	global clicks
	global patch_clicked
	global dic_limits
	global LAST
	figure = plt.gcf()
	if isinstance(event.artist, mt.Text):
		try:
			for t in txt:
				t.remove()
			txt = []
		except:
			pass
	else:
#
# 		try:
# 			txt.remove()
# 		except:
# 			pass
		patch = event.artist
		AXE = patch.axes.get_position()
		dd, axe = [(j,i) for j,i in enumerate(ax) if AXE.p0[1] == i.get_position().p0[1]][0]
		axe.set_zorder(10)
		dd += index_axes
		ice_name = NAME[dd]
		str_ice = [na for na in strains if na in ice_name][0]
		key = [j[0] for j in dic_data[dd].iteritems()
					   if j[1]==[i for i in dic_data[dd].values()
									  if i[2]==patch.get_x()][0]][0]
		texte = "%s\n%s\n%s\n%s\n%s\n%s" %(key,
										   dic_data[dd][key][0][:25],
										   dic_data[dd][key][0][25:50],
										   dic_data[dd][key][8],
										   dic_data[dd][key][9],
										   dic_data[dd][key][7])

		if event.mouseevent.button == 1:

			txt.append(axe.annotate(texte,
									xy=(patch.get_x()+(patch.get_width())/2.,
									dic_data[dd][key][6]), arrowprops = dict(arrowstyle="->"),
									xytext = (patch.get_x(),
											  0.6*len(STRAIN[dd])*dic_data[dd][key][6] if dic_data[dd][key][6]==1 else 1.2*len(STRAIN[dd])*dic_data[dd][key][6]),
									bbox = dict(facecolor='white',  edgecolor='k', alpha=0.85)))

		elif event.mouseevent.button == 3:
			# Every 2 right-clicks, the last 2 are stored as the limits
			# They must be in the same ICE

			if event.mouseevent.key=="shift":

				try:
					for p in dic_patches[ice_name]:
						p.set_edgecolor("k")
						p.set_linewidth(1)
						clicks = []
				except KeyError:
					print "pass"
					pass
				dic_limits = dic_limits[dic_limits.ICE_ID!=ice_name]
				if LAST!=0:
					erase_last_print(LAST)
				print dic_limits
				LAST = len(dic_limits)+1

				del dic_patches[ice_name]

			else:
				clicks.append(key)
				try:
					dic_patches[ice_name].append(patch)
				except KeyError:
					dic_patches[ice_name] = []
					dic_patches[ice_name].append(patch)
				patch.set_edgecolor("red")
				patch.set_linewidth(3)
				if len(clicks)==2:
					first, last = min(clicks), max(clicks)
					lim = pd.DataFrame([[ice_name, first, last]],
					 				    columns=["ICE_ID", "first_gene", "last_gene"])
					dic_limits = pd.concat([dic_limits, lim])
					if LAST!=0:
						erase_last_print(LAST)
					print dic_limits
					LAST = len(dic_limits)+1
					#dic_patches[ice_name] = patch_clicked
					clicks = []




		else:
			pass
	figure.canvas.draw()




############### ---------- Figure ---------- ###############
strains = np.sort(uniq(sum([S.keys() for S in STRAIN],[])))
cm = plt.get_cmap("jet")
cgen = (cm(0.8*i/(len(strains)+1+1)) for i in range(len(strains)+1+3))
next(cgen)
next(cgen)
next(cgen)
col_strain = {}
for s in strains:
	col_strain[s] = next(cgen)

def plot_fig(dic_data,NAME, GC, STRAIN):
	global ax, fig
	fig, ax = plt.subplots(len(dic_data), 1, figsize = (22, 14))
	if isinstance(ax,np.ndarray) == False:
		ax = [ax]
	fig.subplots_adjust(hspace = 0.35)
	#fig.tight_layout(pad = 2.5) #pb Backend
	fig.subplots_adjust(bottom = 0.025, left = 0.035, top = 0.95, right = 0.8)

	ax_barh=[]
	for i in range(len(dic_data)):
		ax_barh.append(plt.axes([ax[i].get_position().get_points()[1][0]+0.01,
								 ax[i].get_position().get_points()[0][1],
								 0.075,
								 ax[i].get_position().get_points()[1][1] - ax[i].get_position().get_points()[0][1]]))
		ax_barh[i].tick_params(top = 'off', left="off", right = 'off', bottom = "off")
		ax_barh[i].set_xticks([])
		ax_barh[i].set_yticks([])
	ax_barh[0].set_title(r"$ \frac{\# gene\ in\ locus}{largest\ locus} $", fontsize=20, va="bottom")



	tRNA = re.compile(r'[\w]*tRNA[\w]*')
	g = []
	pan = []
	bar_rna = []

	for j, dic_d in enumerate(dic_data):

		strain_ice = [na for na in strains if na in NAME[j]][0]

		key = np.sort(dic_d.keys())

		width = np.array([dic_d[k][7] for k in key]) # taille du gène en nt
		position = np.array([dic_d[k][2] for k in key]) # position sur le locus en nt
		heigth = np.array([dic_d[k][6] for k in key]) # 1 ou -1 selon le sens
		hatch = ["///" if dic_d[k][8] != 0 else "xxx" if dic_d[k][9] != 0 else None for k in key ]
		ratiomax = np.max(STRAIN[j].values())
		ratio = np.array([int(STRAIN[j][r]) for r in np.sort(STRAIN[j].keys())])/float(ratiomax)*100.

		g.append(ax[j].bar(position,
		 				   heigth,
						   width,
						   color=col_strain[strain_ice],
						   picker=True))
		for i,gene in enumerate(g[j]):
			gene.set_hatch(hatch[i])
			if dic_d[key[i]][11] != 0 : #RNA

				bar_rna.append(ax[j].bar(dic_d[key[i]][11][2],
								 2*dic_d[key[i]][11][4]*dic_d[key[i]][5],
								 width=dic_d[key[i]][11][3] - dic_d[key[i]][11][2],
								 bottom=dic_d[key[i]][11][4]*(len(STRAIN[j])-1)*dic_d[key[i]][5],
								 color=col_strain[strain_ice],
								 ec=col_strain[strain_ice]))
			if dic_d[key[i]][10] != 0 : #pan
				bot = dic_d[key[i]][6]

				for pan_gene in np.sort(dic_d[key[i]][10]):

					strain_pan = [pg for pg in strains if pg in pan_gene][0]
					pan.append(ax[j].bar(dic_d[key[i]][2],
					 					 dic_d[key[i]][6],
										 dic_d[key[i]][7],
										 color=col_strain[strain_pan],
										 bottom=bot,
										 picker=True))
					bot += dic_d[key[i]][6]
		if dic_d.values()[0][5] == -1:
			ax[j].invert_xaxis()

		#### GC ####

		ax[j].plot(np.linspace(dic_d[key[0]][2],dic_d[key[-1]][3],len(GC[j])), ((np.array(GC[j])-np.min(GC[j])) / np.max(GC[j]-np.min(GC[j]))) * 2*len(STRAIN[j])-len(STRAIN[j]), alpha=0.5, lw=1.75, color = "k",ls="-")

		#### ratio taille locus ####
		ratio_bar = []
		for i,ra in enumerate(ratio):
			ratio_bar.append(ax_barh[j].barh(i, ra, height = 1, color = col_strain[np.sort(STRAIN[j].keys())[i]]))
		ax_barh[j].text(0.2, 0.5, "largest : "+str(int(ratiomax)), transform=ax_barh[j].transAxes)


	for i, a in enumerate(ax):
		for label in a.get_xticklabels():  # make the xtick labels pickable
			label.set_picker(True)

		#a.set_title(NAME[i])
		a.set_ybound( -len(STRAIN[i]) - 0.5, len(STRAIN[i]) + 0.5)
		a.tick_params(top = 'off', left="off", right = 'off')
		a.set_yticks([-len(STRAIN[i]),0,len(STRAIN[i])])
		a.set_yticklabels(["$%.2f$" %(np.min(GC[i])),"$%.2f$" %((np.max(GC[i])+np.min(GC[i]))/2.), "$%.2f$" %np.max(GC[i])])
		a.set_ylabel("%GC")
		#a.set_yticklabels(a.get_yticklabels(), visible = False)

	########### Check boxes ##############


	ax_checkbox = plt.axes([ax[0].get_position().get_points()[1][0] + 0.09, ax[-1].get_position().get_points()[0][1], 0.10,
	ax[0].get_position().get_points()[1][1] - ax[-1].get_position().get_points()[0][1]])

	global check
	check = CheckButtons(ax_checkbox, NAME, [0]*len(NAME))

	for i, (cl, cr, lab) in enumerate(zip(check.lines, check.rectangles, check.labels)) :
		cr.set_y(ax[i].get_position().get_points()[0][1] / ax[0].get_position().get_points()[1][1])
		lab.set_y(cr.get_y() + 0.05)
		cr.set_width(0.15)
		cl[1].set_data([cr.get_x(), cr.get_x() + cr.get_width()], [cr.get_y() + cr.get_height(), cr.get_y()])
		cl[0].set_data([cr.get_x(), cr.get_x() + cr.get_width()], [cr.get_y(), cr.get_y() + cr.get_height()])
		lab.set_size('x-small')
		lab.set_x(cr.get_x() + cr.get_width()+0.01)

	########### Ok button ##############


	ax_ok_button = plt.axes([ax_checkbox.get_position().get_points()[0][0],
					   	 	 ax_checkbox.get_position().get_points()[1][1] + 0.01,
							 0.10,
							 0.025 ])  #left #bottom #width #height
	global ok_button
	ok_button = Button(ax_ok_button, 'OK')
	ok_button.on_clicked(ok_button_callback)

	fig.canvas.mpl_connect('pick_event', onpick1)
	plt.show()

global num_graph
num_graph = 4
global index_axes
global dic_limits
dic_limits = pd.DataFrame()
global dic_patches
dic_patches = {}
if len(dic_data) <= num_graph:
	index_axes = 0
	plot_fig(dic_data,NAME,GC,STRAIN)
else:
	for index_axes in xrange(0,len(dic_data),num_graph):
		plot_fig(dic_data[index_axes:index_axes+num_graph],NAME[index_axes:index_axes+num_graph],GC[index_axes:index_axes+num_graph], STRAIN[index_axes:index_axes+num_graph])

if len(dic_limits>0):
	dic_limits.to_csv("Limits_{}.txt".format(id_locus), sep="\t", index=False)

