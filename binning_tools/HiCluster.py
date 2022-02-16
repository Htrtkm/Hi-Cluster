import sys
import os
import re
from scipy import stats
from math import sqrt
import networkx as nx
from igraph import *
import numpy as np
import matplotlib.pyplot as plt
import copy

def make_scdict(scaffold_fa):
	scdict = {}
	sclength_dict = {}
	length_1000 = 0
	with open(scaffold_fa,"r")as data:
		temporary_seq = ""
		temporary_seq_without_n = ""
		sc_num = 0
		for line in data:
			line_trim = line.strip("\n")
			if line[0] != ">":
				temporary_seq += line
				temporary_seq_without_n += line_trim
			if line[0] == ">":
				scaffold_num = re.search("[0-9]{1,8}",line_trim)
				scaffold_num_int = int(scaffold_num.group())
				if temporary_seq == "":
					sc_num = scaffold_num_int
				if temporary_seq != "":
					scdict.setdefault(sc_num,temporary_seq)
					sclength_dict.setdefault(sc_num,len(temporary_seq_without_n))
					sc_num = scaffold_num_int
					if len(temporary_seq_without_n) >= 1000:
						length_1000 += len(temporary_seq_without_n)
				temporary_seq = ""
				temporary_seq_without_n = ""
		if temporary_seq != "":
			scdict.setdefault(sc_num,temporary_seq)
			sclength_dict.setdefault(sc_num,len(temporary_seq_without_n))
		temporary_seq = ""
		temporary_seq_without_n = ""
	print(length_1000)
	return scdict,sclength_dict

def read_gml(gmlfile):
	id_label = {}
	edgelist = {}
	temporary_edgelist = []
	with open(gmlfile,"r")as data:
		node_id = 0
		node_label = 0
		edge_source = 0
		edge_target = 0
		edge_weight = 0.0
		for line in data:
			line_trim = line.strip("\n")
			table = line_trim.split()
			if "id" in line:
				node_id = int(table[1])
			if "label" in line:
				node_label = int(table[1].strip("\""))
				id_label.setdefault(node_id,node_label)
			if "source" in line:
				edge_source = int(table[1])
			if "target" in line:
				edge_target = int(table[1])
			if "weight" in line:
				edge_weight = float(table[1])
				temporary_edgelist.append([edge_source,edge_target,edge_weight])
	temporary_edgelist.sort(key=lambda x:-x[2])
	for i in range(len(temporary_edgelist)):
		edgelist.setdefault(i,temporary_edgelist[i])
	return id_label,edgelist

def select_edges(edgelist,recruit_rate):
	id_edge = {}
	id_edge_selected ={}
	remove_edges = {}
	for k,v in edgelist.items():
		if v[0] in id_edge:
			id_edge[v[0]].append(k)
		id_edge.setdefault(v[0],[k])
		if v[1] in id_edge:
			id_edge[v[1]].append(k)
		id_edge.setdefault(v[1],[k])
	rate_recruited_edge = []
	id_edge_sorted = sorted(id_edge.items(),key=lambda x:x[0])
	id_edge_dict_sorted = {}
	for i in range(len(id_edge_sorted)):
		id_edge_dict_sorted.setdefault(id_edge_sorted[i][0],id_edge_sorted[i][1])
	for k,v in id_edge_dict_sorted.items():
		weights_recruited = 0.0
		weights_all = 0.0
		switch = 0
		removes = 0
		for i in range(len(v)):
			weights_all += edgelist[v[i]][2]
		for i in range(len(v)):
			if switch == 0:
				weights_recruited += edgelist[v[i]][2]
				rate_recruited_edge.append(v[i])
				if weights_recruited/weights_all >= float(recruit_rate):
					switch = 1
					continue
			if switch == 1:
				removes += 1
		if (len(v)-removes)/len(v) <= 0.4:
			v = rate_recruited_edge
			id_edge_selected.setdefault(k,rate_recruited_edge)
		if (len(v)-removes)/len(v) > 0.4:
			id_edge_selected.setdefault(k,v)
		rate_recruited_edge = []
	temporary_edgelist = {}
	temporary_edgelist2 = {}
	temporary_edgelist3 = {}
	edgelist_recruited = {}
	for k,v in id_edge_selected.items():
		for i in range(len(v)):
			if v[i] in temporary_edgelist:
				temporary_edgelist[v[i]] += 1
			temporary_edgelist.setdefault(v[i],1)
	for k,v in temporary_edgelist.items():
		if v == 2:
			temporary_edgelist2.setdefault(k,0)
	temporary_edgelist3 = sorted(temporary_edgelist2.items(),key=lambda x:x[0])
	for i in range(len(temporary_edgelist3)):
		edgelist_recruited.setdefault(temporary_edgelist3[i][0],0)
	return id_edge_selected,edgelist_recruited

def make_gml(id_edge_selected,id_label,edgelist,edgelist_recruited,name):
	first_gml = "first_" + name + ".gml"
	with open(first_gml,"w")as f:
		nodes = ""
		edges = ""
		id2_id = {}
		id_id2 = {}
		id2 = 0
		idlists = []
		for k,v in id_edge_selected.items():
			idlists.append(k)
		idlists.sort()
		for i in range(len(idlists)):
			id2_id.setdefault(id2,idlists[i])
			id_id2.setdefault(idlists[i],id2)
			nodes += "\tnode\n\t[\n\t\tid " + str(id2) + "\n\t\tlabel \"" + str(id_label[idlists[i]]) + "\"\n\t]\n"
			id2 += 1
		for k,v in edgelist_recruited.items():
			edges += "\tedge\n\t[\n\t\tsource " + str(id_id2[edgelist[k][0]]) + "\n\t\ttarget " + str(id_id2[edgelist[k][1]]) + "\n\t\tweight " + str(edgelist[k][2]) + "\n\t]\n"
		f.write("Creator \"Takumi Hattori\"" + "\n" + "graph" + "\n" + "[" + "\n")
		f.write(nodes)
		f.write(edges)
		f.write("]")
	return id2_id

def select_bins(binningresult1,sclength_dict,id_label,idfirst_id,scdict):
	normalbins = {}
	hugebins = {}
	smallbins = {}
	normalbin_length = {}
	hugebin_length = {}
	id_normallength = {}
	id_hugelength = {}
	bin_length_sum = 0
	for i in range(len(binningresult1)):
		bin_length = 0
		contigs_in_the_bin = []
		for j in range(len(binningresult1[i])):
			bin_length += sclength_dict[id_label[idfirst_id[int(binningresult1[i][j])]]]
			contigs_in_the_bin.append(idfirst_id[int(binningresult1[i][j])])
		if bin_length < 8000000 :
			binname = "bin-"  + str(i) + "-normal.fa"
			if bin_length >= 1000000:
				print(binname + " " + str(bin_length))
			normalbins.setdefault(binname,contigs_in_the_bin)
			normalbin_length.setdefault(binname,bin_length)
			#if bin_length >= 500000:
				#with open("contigsum.txt","w")as f:
					#f.write(binname + "\n")
					#f.write("bin_length :" + str(bin_length) + "\n")
					#f.write("contig sum :" + str(len(contigs_in_the_bin)) + "\n\n")
		if bin_length >= 8000000:
			binname = "bin-" + str(i) + "-huge.fa"
			print(binname + " " + str(bin_length))
			hugebins.setdefault(binname,contigs_in_the_bin)
			hugebin_length.setdefault(binname,bin_length)
		bin_length_sum += bin_length
	for k,v in normalbins.items():
		for i in range(len(v)):
			id_normallength.setdefault(v[i],normalbin_length[k])
	for k,v in hugebins.items():
		for i in range(len(v)):
			id_hugelength.setdefault(v[i],hugebin_length[k])
	print("first length :" + str(bin_length_sum))
	return normalbins,hugebins,id_normallength,id_hugelength

def make_contig_firstbins(normalbins_first,hugebins_first):
	id_firstbin = {}
	for k,v in normalbins_first.items():
		for i in range(len(v)):
			id_firstbin.setdefault(v[i],k)
	for k,v in hugebins_first.items():
		for i in range(len(v)):
			id_firstbin.setdefault(v[i],k)
	return id_firstbin

def makegml_normal(normalbins_first,edgelist_recruited,edgelist,id_label,scdict,sclength_dict,id_normallength,name):
	edgelist_normal = {}
	id_normalbin = {}
	id_normalbin_list = []
	id_edge = {}
	removes = {}
	for k,v in normalbins_first.items():
		for i in range(len(v)):
			id_normalbin_list.append(v[i])
	id_normalbin_list.sort()
	for i in range(len(id_normalbin_list)):
		id_normalbin.setdefault(id_normalbin_list[i],0)
	for k,v in edgelist_recruited.items():
		if edgelist[k][0] in id_normalbin and edgelist[k][1] in id_normalbin:
			edgelist_normal.setdefault(k,edgelist[k])
			if edgelist[k][0] in id_edge:
				id_edge[edgelist[k][0]].append(k)
			id_edge.setdefault(edgelist[k][0],[k])
			if edgelist[k][1] in id_edge:
				id_edge[edgelist[k][1]].append(k)
			id_edge.setdefault(edgelist[k][1],[k])
	for k,v in id_edge.items():
		sum_weight = 0.0
		sum_recruit_weight = 0.0
		switch = 0
		recruit_rate = 1.0
		if id_normallength[k] >= 2000000 and id_normallength[k] < 2500000:
			recruit_rate = 0.99
		if id_normallength[k] >= 2500000 and id_normallength[k] < 3000000:
			recruit_rate = 0.98
		if id_normallength[k] >= 3000000 and id_normallength[k] < 3500000:
			recruit_rate = 0.96
		if id_normallength[k] >= 3500000 and id_normallength[k] < 4000000:
			recruit_rate = 0.96
		if id_normallength[k] >= 4000000 and id_normallength[k] < 5000000:
			recruit_rate = 0.92
		if id_normallength[k] >= 5000000:
			recruit_rate = 0.9
		for i in range(len(v)):
			sum_weight += edgelist[v[i]][2]
		for i in range(len(v)):
			sum_recruit_weight += edgelist[v[i]][2]
			if switch == 0:
				if sum_recruit_weight/sum_weight >= recruit_rate:
					switch = 1
					continue
			if switch == 1:
				removes.setdefault(v[i],0)
	nodes = ""
	edges = ""
	idnormal_id = {}
	id_idnormal = {}
	edges_normal = {}
	idnormal = 0
	normal_gml = "normal_" + name + ".gml"
	with open(normal_gml,"w")as f:
		for k,v in id_normalbin.items():
			idnormal_id.setdefault(idnormal,k)
			id_idnormal.setdefault(k,idnormal)
			nodes += "\tnode\n\t[\n\t\tid " + str(idnormal) + "\n\t\tlabel \"" + str(id_label[k]) + "\"\n\t]\n"
			idnormal += 1
		for k,v in edgelist_normal.items():
			if k not in removes:
				edges_normal.setdefault(k,v)
				edges += "\tedge\n\t[\n\t\tsource " + str(id_idnormal[v[0]]) + "\n\t\ttarget " + str(id_idnormal[v[1]]) + "\n\t\tweight " + str(v[2]) + "\n\t]\n"
		f.write("Creator \"Takumi Hattori\"" + "\n" + "graph" + "\n" + "[" + "\n")
		f.write(nodes)
		f.write(edges)
		f.write("]")
	gn = Graph.Read_GML(normal_gml)
	binningresult_normal = gn.community_infomap()
	smallbins_normal = {}
	normalbins_second = {}
	rescuecontig = {}
	for i in range(len(binningresult_normal)):
		bin_length = 0
		contigs_in_the_bin = []
		for j in range(len(binningresult_normal[i])):
			bin_length += sclength_dict[id_label[idnormal_id[int(binningresult_normal[i][j])]]]
			contigs_in_the_bin.append(idnormal_id[int(binningresult_normal[i][j])])
		binname = str(i) + "_normal.fa"
		with open(binname,"w")as f:
			for j in range(len(binningresult_normal[i])):
			 	f.write(">" + str(id_label[idnormal_id[int(binningresult_normal[i][j])]]) + "\n")
			 	f.write(scdict[id_label[idnormal_id[int(binningresult_normal[i][j])]]])
		normalbins_second.setdefault(binname,contigs_in_the_bin)
		if bin_length <= 2000000:
			smallbins_normal.setdefault(binname,contigs_in_the_bin)
		for i in range(len(contigs_in_the_bin)):
			if sclength_dict[id_label[contigs_in_the_bin[i]]] >= 10000000:
				#rescue contig wo sakujo
				rescuecontig.setdefault(contigs_in_the_bin[i],0)
	os.remove(normal_gml)
	return normalbins_second,smallbins_normal,edges_normal,rescuecontig

def repeat_localbinning(hugebins_first,edgelist_recruited,edgelist,sclength_dict,id_label,name):
	smallbins_huge_all = {}
	huge_to_normal_all ={}
	edges_huge_all ={}
	hugebin_num = 0
	for k,v in hugebins_first.items():
		try_sum = 0
		hugebins_info = [[k,v]]
		print(str(k))
		while hugebins_info != []:
			huge_bin_info = local_binning(edgelist_recruited,edgelist,hugebins_info[0][0],hugebins_info[0][1],sclength_dict,id_label,hugebin_num,name)
			smallbins_huge = huge_bin_info[0]
			huge_to_normal = huge_bin_info[1]
			still_big_one = huge_bin_info[2]
			edges_huge = huge_bin_info[3]
			the_num = huge_bin_info[4]
			hugebin_num += the_num
			contigs_in_still_big_one = {}
			if still_big_one != {}:
				for l,m in still_big_one.items():
					for i in range(len(m)):
						contigs_in_still_big_one.setdefault(m[i],l)
			for l,m in smallbins_huge.items():
				smallbins_huge_all.setdefault(l,m)
			for l,m in huge_to_normal.items():
				huge_to_normal_all.setdefault(l,m)
			for l,m in edges_huge.items():
				if m[0] not in contigs_in_still_big_one and m[1] not in contigs_in_still_big_one:
					edges_huge_all.setdefault(l,m)
			print(hugebins_info[0][0])
			hugebins_info.pop(0)
			for l,m in still_big_one.items():
				hugebins_info.append([l,m])
			try_sum += 1
			if try_sum >= 10:
				for l,m in huge_to_normal.items():
					huge_to_normal_all.setdefault(hugebins_info[0][0],hugebins_info[0][1])
				for l,m in edges_huge.items():
					if m[0] in contigs_in_still_big_one and m[1] in contigs_in_still_big_one:
						edges_huge_all.setdefault(l,m)
				hugebins_info = []
	return smallbins_huge_all,huge_to_normal_all,edges_huge_all

def merge_edges(edges_normal,edges_huge_all,smallbins_normal,smallbins_huge_all):
	edges_used = {}
	for k,v in edges_normal.items():
		edges_used.setdefault(k,v)
	for k,v in edges_huge_all.items():
		edges_used.setdefault(k,v)
	smallbins_merged = {}
	for k,v in smallbins_normal.items():
		smallbins_merged.setdefault(k,v)
	for k,v in smallbins_huge_all.items():
		smallbins_merged.setdefault(k,v)
	return edges_used,smallbins_merged

def merge_smallbins(edges_used,smallbins_merged,name):
	id_bin = {}
	bin_weights = {}
	for k,v in smallbins_merged.items():
		for i in range(len(v)):
			id_bin.setdefault(v[i],k)
	for k,v in smallbins_merged.items():
		bin_weights.setdefault(k,{})	
	for k,v in edges_used.items():
		if v[0] in id_bin and v[1] in id_bin:
			if id_bin[v[0]] == id_bin[v[1]]:
				if id_bin[v[0]] in bin_weights[id_bin[v[0]]]:
					bin_weights[id_bin[v[0]]][id_bin[v[0]]] += v[2]
				bin_weights[id_bin[v[0]]].setdefault(id_bin[v[0]],v[2])
			if id_bin[v[0]] != id_bin[v[1]]:
				if id_bin[v[1]] in bin_weights[id_bin[v[0]]]:
					bin_weights[id_bin[v[0]]][id_bin[v[1]]] += v[2]
				bin_weights[id_bin[v[0]]].setdefault(id_bin[v[1]],v[2])
				if id_bin[v[0]] in bin_weights[id_bin[v[1]]]:
					bin_weights[id_bin[v[1]]][id_bin[v[0]]] += v[2]
				bin_weights[id_bin[v[1]]].setdefault(id_bin[v[0]],v[2])
	binslist = {}
	for k,v in bin_weights.items():
		weight_sum = 0.0
		inner_weight = 0.0
		recruit_weight = 0.0
		v_sorted = sorted(v.items(),key=lambda x:-x[1])
		for i in range(len(v_sorted)):
			weight_sum += v_sorted[i][1]
			if k == v_sorted[i][0]:
				inner_weight = v_sorted[i][1]
		if weight_sum == 0:
			continue
		if (weight_sum - inner_weight)/weight_sum <= 0.15:
			continue
		if k in smallbins_merged:
			binslist.setdefault(k,[])
		if k not in smallbins_merged:
			continue
		recruit_sum = 0
		for i in range(len(v_sorted)):
			if k != v_sorted[i][0]:
				if v_sorted[i][0] in smallbins_merged :
					binslist[k].append([v_sorted[i][0],v_sorted[i][1]])
					recruit_weight += v_sorted[i][1]
					recruit_sum += 1
					if recruit_sum == 2:
						break
	small_gml = "small_" + name + ".gml"
	with open(small_gml,"w")as f:
		f.write("Creator \"Takumi Hattori\"\ngraph\n[\n")
		id_bins = {}
		bins_id = {}
		binid = 0
		bins_recruited = {}
		edges_recruited = {}
		edges_recruited_huge = {}
		for k,v in binslist.items():
			for i in range(len(v)):
				if v[i][0] in binslist:
					bins_recruited.setdefault(v[i][0],0)
					bins_recruited.setdefault(k,0)
					bin_num1 = int(re.search("[0-9]{1,8}",k).group())
					bin_num2 = int(re.search("[0-9]{1,8}",v[i][0]).group())
					if "normal" in k :
						if bin_num1 >= bin_num2:
							if v[i][0] + "," + k in edges_recruited:
								edges_recruited[v[i][0] + "," + k] += 1
							edges_recruited.setdefault(v[i][0] + "," + k,1)
						if bin_num2 > bin_num1:
							if k + "," + v[i][0] in edges_recruited:
								edges_recruited[k + "," + v[i][0]] += 1
							edges_recruited.setdefault(k + "," + v[i][0],1)
					if "huge" in k :
						if bin_num1 >= bin_num2:
							if v[i][0] + "," + k in edges_recruited_huge:
								edges_recruited_huge[v[i][0] + "," + k] += 1
							edges_recruited_huge.setdefault(v[i][0] + "," + k,1)
						if bin_num2 > bin_num1:
							if k + "," + v[i][0] in edges_recruited_huge:
								edges_recruited_huge[k + "," + v[i][0]] += 1
							edges_recruited_huge.setdefault(k + "," + v[i][0],1)
		for k,v in bins_recruited.items():
			f.write("\tnode\n\t[\n\t\tid " + str(binid) + "\n\t\tlabel \"" + k + "\"\n\t]\n")
			id_bins.setdefault(binid,k)
			bins_id.setdefault(k,binid)
			binid += 1
		for k,v in edges_recruited.items():
			if v == 2 or v == 1:
				table = k.split(",")
				source = table[0]
				target = table[1]
				f.write("\tedge\n\t[\n\t\tsource " + str(bins_id[source]) + "\n\t\ttarget " + str(bins_id[target]) + "\n\t]\n")
		for k,v in edges_recruited_huge.items():
			if v == 2 or v == 1:
				table = k.split(",")
				source = table[0]
				target = table[1]
				f.write("\tedge\n\t[\n\t\tsource " + str(bins_id[source]) + "\n\t\ttarget " + str(bins_id[target]) + "\n\t]\n")
		f.write("]")
	gs = Graph.Read_GML(small_gml)
	os.remove(small_gml)
	binningresult_small = gs.community_infomap()
	for i in range(len(binningresult_small)):
		merged_bins = []
		for j in range(len(binningresult_small[i])):
			merged_bins.append(id_bins[binningresult_small[i][j]])
		#print(merged_bins)
	binnum_small = 0
	mergedbins = {}
	for i in range(len(binningresult_small)):
		binname = str(binnum_small) + "-small.fa"
		binnum_small += 1
		if len(binningresult_small[i]) >= 2:
			with open(binname,"w")as f:
				for j in range(len(binningresult_small[i])):
					mergedbins.setdefault(id_bins[binningresult_small[i][j]],binname)
					bin_name = id_bins[binningresult_small[i][j]]
					with open(bin_name,"r")as data:
						for line in data:
							f.write(line)
					os.remove(id_bins[binningresult_small[i][j]])
	return mergedbins

def make_bin_contigs(normalbins_second,huge_to_normal_all,mergedbins):
	bin_contigs = {}
	for k,v in normalbins_second.items():
		if k not in mergedbins:
			bin_contigs.setdefault(k,v)
		if k in mergedbins:
			if mergedbins[k] in bin_contigs:
				bin_contigs[mergedbins[k]].extend(v)
			bin_contigs.setdefault(mergedbins[k],v)
	for k,v in huge_to_normal_all.items():
		if k not in mergedbins:
			bin_contigs.setdefault(k,v)
		if k in mergedbins:
			if mergedbins[k] in bin_contigs:
				bin_contigs[mergedbins[k]].extend(v)
			bin_contigs.setdefault(mergedbins[k],v)
	return bin_contigs

def third_binning(edges_used,sclength_dict,id_label,bin_contigs):
	clear_edges = {}
	id_edges_all = {}
	for k,v in edges_used.items():
		if v[0] in id_edges_all:
			id_edges_all[v[0]].append(k)
		id_edges_all.setdefault(v[0],[k])
		if v[1] in id_edges_all:
			id_edges_all[v[1]].append(k)
		id_edges_all.setdefault(v[1],[k])
	for k,v in bin_contigs.items():
		bin_length = 0
		for i in range(len(v)):
			bin_length += sclength_dict[id_label[v[i]]]
		if bin_length >= 5000000 or len(v) >= 10:
			temporary_edges = {}
			id_edges = {}
			for l,m in edges_used.items():
				if m[0] in v:
					temporary_edges.setdefault(l,m)
					if m[0] in id_edges:
						id_edges[m[0]].append(l)
					id_edges.setdefault(m[0],[l])
				if m[1] in v:
					temporary_edges.setdefault(l,m)
					if m[1] in id_edges:
						id_edges[m[1]].append(l)
					id_edges.setdefault(m[1],[l])
			edge_score = {}
			for l,m in id_edges.items():
				if sclength_dict[id_label[l]] >= 800000:
					edge_candidates = {}
					temporary_edge_score = {}
					neighborhood = {}
					for j in range(len(m)):
						if temporary_edges[m[j]][0] == l:
							edge_candidates.setdefault(temporary_edges[m[j]][1],[m[j],0])
							neighborhood.setdefault(temporary_edges[m[j]][1],id_edges_all[temporary_edges[m[j]][1]])
						if temporary_edges[m[j]][1] == l:
							edge_candidates.setdefault(temporary_edges[m[j]][0],[m[j],0])
							neighborhood.setdefault(temporary_edges[m[j]][0],id_edges_all[temporary_edges[m[j]][0]])
					for k2,v2 in neighborhood.items():
						for j in range(len(v2)):
							if v2[j] in temporary_edges:
								if temporary_edges[v2[j]][0] == k2:
									if temporary_edges[v2[j]][1] not in neighborhood:
										if temporary_edges[v2[j]][1] != l:
											edge_candidates[k2][1] += 1
								if temporary_edges[v2[j]][1] == k2:
									if temporary_edges[v2[j]][0] not in neighborhood:
										if temporary_edges[v2[j]][0] != l:
											edge_candidates[k2][1] += 1
					for l2,m2 in edge_candidates.items():
						temporary_edge_score.setdefault(m2[0],m2[1])
					for l2,m2 in temporary_edge_score.items():
						edge_score.setdefault(l2,0)
						if m2 == 0:
							edge_score[l2] += 1
			for l,m in temporary_edges.items():
				if l not in edge_score:
					clear_edges.setdefault(l,m)
				if l in edge_score:
					if edge_score[l] >= 1:
						clear_edges.setdefault(l,m)
		else:
			temporary_edges = {}
			for l,m in edges_used.items():
				if m[0] in v:
					temporary_edges.setdefault(l,m)
				if m[1] in v:
					temporary_edges.setdefault(l,m)
			for l,m in temporary_edges.items():
				clear_edges.setdefault(l,m)
	print(len(edges_used))
	print(len(edges_used)-len(clear_edges))
	return clear_edges

def third_binning_latter(clear_edges,bin_contigs):
	id_id3 = {}
	id3_id = {}
	with open("third.gml","w")as f:
		nodes_info = ""
		edges = ""
		id3 = 0
		for k,v in bin_contigs.items():
			for i in range(len(v)):
				id3_id.setdefault(id3,v[i])
				id_id3.setdefault(v[i],id3)
				nodes_info += "\tnode\n\t[\n\t\tid " + str(id3) + "\n\t\tlabel \"" + str(id_label[id3_id[id3]]) + "\"\n\t]\n"
				id3 += 1
		for k,v in clear_edges.items():
			if v[0] in id_id3 and v[1] in id_id3:
				edges += "\tedge\n\t[\n\t\tsource " + str(id_id3[v[0]]) + "\n\t\ttarget " + str(id_id3[v[1]]) + "\n\t\tweight " + str(v[2]) + "\n\t]\n"
		f.write("Creator \"Takumi Hattori\"\ngraph\n[\n")
		f.write(nodes_info)
		f.write(edges)
		f.write("]")
	g3 = Graph.Read_GML("third.gml")
	os.remove("third.gml")
	result3 = g3.community_infomap()
	bin3_contigs = {}
	for i in range(len(result3)):
		name3 = str(i) + "_thrid.fa"
		bin_length = 0
		contigs_in_the_bin = []
		for j in range(len(result3[i])):
			bin_length += sclength_dict[id_label[id3_id[int(result3[i][j])]]]
			contigs_in_the_bin.append(id3_id[int(result3[i][j])])
		if bin_length >= 100000:
			bin3_contigs.setdefault(name3,contigs_in_the_bin)
	return bin3_contigs

def make_third_bins(bin3_contigs,name,sclength_dict,id_label,scdict):
	third_dir = "trimedbins_" + name
	os.mkdir(third_dir)
	os.chdir(third_dir)
	sum_length = 0
	for k,v in bin3_contigs.items():
		bin_length = 0
		for i in range(len(v)):
			bin_length += sclength_dict[id_label[v[i]]]
		sum_length += bin_length
		if bin_length >= 500000:
			with open(k,"w")as f:
				for i in range(len(v)):
					f.write(">" + str(id_label[v[i]]) + "\n")
					f.write(scdict[id_label[v[i]]])
	print(sum_length)

def local_binning(edgelist_recruited,edgelist,binname_original,contigs_in_original,sclength_dict,id_label,hugebin_num,name):
	contigs = {}
	contigs_list = []
	id_edge = {}
	edges_in_original = {}
	removes = {}
	for i in range(len(contigs_in_original)):
		contigs_list.append(contigs_in_original[i])
	contigs_list.sort()
	for i in range(len(contigs_list)):
		contigs.setdefault(contigs_list[i],0)
	for k,v in edgelist_recruited.items():
		if edgelist[k][0] in contigs and edgelist[k][1] in contigs:
			edges_in_original.setdefault(k,edgelist[k])
			if edgelist[k][0] in id_edge:
				id_edge[edgelist[k][0]].append(k)
			id_edge.setdefault(edgelist[k][0],[k])
			if edgelist[k][1] in id_edge:
				id_edge[edgelist[k][1]].append(k)
			id_edge.setdefault(edgelist[k][1],[k])
	for k,v in id_edge.items():
		sum_weight = 0.0
		sum_recruit_weight = 0.0
		switch = 0
		recruit_rate = 1.0
		temporary_removes = {}
		remove_num = 0
		for i in range(len(v)):
			sum_weight += edgelist[v[i]][2]
		for i in range(len(v)):
			sum_recruit_weight += edgelist[v[i]][2]
			if switch == 0:
				if sum_recruit_weight/sum_weight >= 0.9:
					switch = 1
					continue
			if switch == 1:
				remove_num += 1
				temporary_removes.setdefault(v[i],0)
		if (len(v)-remove_num)/len(v) <= 1:
			for l,m in temporary_removes.items():
				removes.setdefault(l,m)
	nodes = ""
	edges = ""
	idhuge_id = {}
	id_idhuge = {}
	edges_huge = {}
	idhuge = 0
	huge_gml = "huge_" + name + ".gml"
	with open(huge_gml,"w")as f:
		for k,v in contigs.items():
			idhuge_id.setdefault(idhuge,k)
			id_idhuge.setdefault(k,idhuge)
			nodes += "\tnode\n\t[\n\t\tid " + str(idhuge) + "\n\t\tlabel \"" + str(id_label[k]) + "\"\n\t]\n"
			idhuge += 1
		for k,v in edges_in_original.items():
			if k not in removes:
				edges_huge.setdefault(k,v)
				edges += "\tedge\n\t[\n\t\tsource " + str(id_idhuge[v[0]]) + "\n\t\ttarget " + str(id_idhuge[v[1]]) + "\n\t\tweight " + str(v[2]) + "\n\t]\n"
		f.write("Creator \"Takumi Hattori\"" + "\n" + "graph" + "\n" + "[" + "\n")
		f.write(nodes)
		f.write(edges)
		f.write("]")
	gh = Graph.Read_GML(huge_gml)
	binningresult_huge = gh.community_infomap()
	smallbins_huge = {}
	huge_to_normal = {}
	still_big_one = {}
	the_number = 0
	for i in range(len(binningresult_huge)):
		bin_length = 0
		the_number += 1
		contigs_in_the_bin = []
		for j in range(len(binningresult_huge[i])):
			bin_length += sclength_dict[id_label[idhuge_id[int(binningresult_huge[i][j])]]]
			contigs_in_the_bin.append(idhuge_id[int(binningresult_huge[i][j])])
		if bin_length >= 100000:
			binname = str(the_number + hugebin_num) + "__huge.fa"
			if bin_length < 8000000:
				with open(binname,"w")as f:
					for j in range(len(binningresult_huge[i])):
						f.write(">" + str(id_label[idhuge_id[int(binningresult_huge[i][j])]]) + "\n")
						f.write(scdict[id_label[idhuge_id[int(binningresult_huge[i][j])]]])
			if bin_length < 8000000:
				huge_to_normal.setdefault(binname,contigs_in_the_bin)
			if bin_length <= 2000000:
				smallbins_huge.setdefault(binname,contigs_in_the_bin)
			if bin_length >= 8000000:
				still_big_one.setdefault(binname,contigs_in_the_bin)
				print("still big one exists")
				print(bin_length)
	os.remove(huge_gml)
	return smallbins_huge,huge_to_normal,still_big_one,edges_huge,the_number

args = sys.argv

#ready for first binning
result_scdict = make_scdict(args[2])
scdict = result_scdict[0]
sclength_dict = result_scdict[1]

result_readgml = read_gml(args[1])
id_label = result_readgml[0]
edgelist = result_readgml[1]

result_select_edges = select_edges(edgelist,args[3])
id_edge_selected = result_select_edges[0]
edgelist_recruited = result_select_edges[1]

idfirst_id = make_gml(id_edge_selected,id_label,edgelist,edgelist_recruited,args[4])


#first_binning
first_gml = "first_" + args[4] + ".gml"
g1 = Graph.Read_GML(first_gml)
binningresult1 = g1.community_infomap()
print("finish first binning")
os.remove(first_gml)

result_select_bins = select_bins(binningresult1,sclength_dict,id_label,idfirst_id,scdict)
normalbins_first = result_select_bins[0]
hugebins_first = result_select_bins[1]
id_normallength = result_select_bins[2]
id_firstbin = make_contig_firstbins(normalbins_first,hugebins_first)

#second_binning
directory = "bins" + args[4]
os.mkdir(directory)
os.chdir(directory)

result_makegml_normal = makegml_normal(normalbins_first,edgelist_recruited,edgelist,id_label,scdict,sclength_dict,id_normallength,args[4])
normalbins_second = result_makegml_normal[0]
smallbins_normal = result_makegml_normal[1]
edges_normal = result_makegml_normal[2]
rescuecontig = result_makegml_normal[3]

result_repeat_localbinning = repeat_localbinning(hugebins_first,edgelist_recruited,edgelist,sclength_dict,id_label,args[4])
smallbins_huge_all = result_repeat_localbinning[0]
huge_to_normal_all = result_repeat_localbinning[1]
edges_huge_all = result_repeat_localbinning[2]
print("finish second binning")


#merge_smallbins
edges_used_smallbins_merged = merge_edges(edges_normal,edges_huge_all,smallbins_normal,smallbins_huge_all)
edges_used = edges_used_smallbins_merged[0]
smallbins_merged = edges_used_smallbins_merged[1]
mergedbins = merge_smallbins(edges_used,smallbins_merged,args[4])

#remove suspectable contigs
os.chdir('../')
bin_contigs = make_bin_contigs(normalbins_second,huge_to_normal_all,mergedbins)
clear_edges = third_binning(edges_used,sclength_dict,id_label,bin_contigs)
bin3_contigs = third_binning_latter(clear_edges,bin_contigs)
make_third_bins(bin3_contigs,str(args[4]),sclength_dict,id_label,scdict)
print("finish third binning")
"""
args1 ->gml file
args2 ->fasta file
args3 ->0.96
args4 ->name
"""
