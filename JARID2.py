#!/usr/bin/python

##############################################################################################
import os
os.environ['HOME'] = "./" #required for write rights for matplotlib in www-user environment	 #	
##############################################################################################

import matplotlib as plt
plt.use("Agg")
import warnings
import cPickle
import sys
import argparse
import pylab
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import collections
import numpy
import h5py
import rpy2.rlike.container as rlc



#Global vars :
#-------------

path_to_use = './'

#-------------


#################################################### Functions ############################################

#We do not want this depreciated warning
def fxn():
	warnings.warn("deprecated", DeprecationWarning)
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	fxn()
	pass
#
def import_sample_data(_file,sep='\t',verbose=False, line_to_skip=0):
	"""docstring for import_sample_data
	import data from tab separated file containing sample description. 
	"""
	f=open(_file,'r')
	t=f.readlines()
	f.close()
	data=[]
	for sample in t[line_to_skip:]:
		data.append( [i.strip() for i in sample.strip().split(sep)] )
		if verbose:
			print '#%s#'% data[-1] 
	return data
	pass
#

def beeswarm(data, names,out='beeswarm.pdf',title='',plot_median=False, use_mean=False, use_log=False):
	"""docstring for beeswarm
	plots a beeswram plot.
	data is an array or values
	names is the corresponding names for each data values.
	
	"""
	import pylab
	if use_mean:
		meam_method = numpy.mean
	else:
		meam_method = numpy.median
		
	in_gene = title
	graphics=importr('graphics')
	beeswarm= importr('beeswarm')
	grdevices = importr('grDevices')
	all_classes_order=[]
	all_classes={}
	for index,classe in enumerate(names):
		if classe not in all_classes:
			all_classes_order.append(classe)
			all_classes[classe]=[]
			all_classes[classe].append(index)
		else:
			all_classes[classe].append(index)
	
	
	nn=robjects.StrVector(numpy.array(names))
	fcdf=robjects.FloatVector(numpy.array(data))
	cool=robjects.DataFrame({'Expression': fcdf , 'Class':nn})
	Classes= robjects.r.factor(cool.rx2(2), levels = robjects.r.unique(cool.rx2(2)), ordered = 1)
	robjects.globalenv['Classes']=Classes
	
	a=beeswarm.beeswarm( robjects.r("Expression ~ Classes") , data = cool, method = "center",pch=1, bg=robjects.r("rainbow(6)") ,col = robjects.r("rainbow(8)") , las =2 )
	opt={'main':in_gene}
	graphics.title( **opt)
	grdevices.dev_off()
	
	x=numpy.array(a[0])
	y=numpy.array(a[1])
	c=numpy.array(a[3])
	
	fig=pylab.figure()
	
	if use_log:
		pylab.yscale('log', basey=2)
		pylab.ylabel('log2 expression')
		
		
	for i,j in enumerate(x): 
		pylab.plot(x[i],numpy.exp2(y[i]),'o',color=c[i][0:-2],alpha=.7) # R adds FF at the end of the color string, which is bad.
		
	if plot_median:
		all_classes_order=[]
		all_classes={}
#		names = numpy.array( [ i.split(')')[0]+')'	for i in a.rownames ] ) # F#? R renames things...
		for index,classe in enumerate(names):
			if classe not in all_classes:
				all_classes_order.append(classe)
				all_classes[classe]=[]
				all_classes[classe].append(index)
			else:
				all_classes[classe].append(index)
		
		for	 i in range(len (all_classes_order)):
			z= numpy.array(all_classes[all_classes_order[i]]) 
			m=meam_method(numpy.exp2(y[z]))
			pylab.plot([i+1-.3, i+1+.3],[m,m],'b',linewidth=1.5, alpha=.7)
#		pylab.plot([1,2],[0,0],'b',linewidth=1.5,alpha=.5)
	pylab.xticks(range(len(all_classes_order)+1), numpy.concatenate([[''],all_classes_order],axis=0), rotation=90,size=9)
	fig.autofmt_xdate(bottom=0.18)	
		
	pylab.title(in_gene)
	pylab.savefig(out)
	pylab.close()



def pickle_load(in_file):
	"""docstring for pickle_load
	load a pickled variable into memory.
	input is the file name and path.
	"""
	x=None
	if os.path.exists(in_file):
		pkl_file = open(in_file, 'rb')
		x=cPickle.load(pkl_file)
		pkl_file.close()
	else:
		print 'File %s could not be found ! ' % (in_file)
	
	return x
	pass

def find_probes(gene,organism):
	"""docstring for find_probes
		loads the annotation DBI R package converted by the generate_annotation_dics.py script.
		loops through all possible annotations, and the fuction hopefully returns the correct probes ;)
	"""
	gene=gene.lower()
	
	if 0: #slow method.
		p=None
		if organism == 'human':
			annot_to_probe = pickle_load(path_to_use + '/beeswarn_data/h_annot_to_probe.pkl')
		else:
			annot_to_probe = pickle_load(path_to_use + '/beeswarn_data/m_annot_to_probe.pkl')
			
			
		for k in annot_to_probe:
			p = annot_to_probe[k].get(gene)
			if p != None:
				if not WEB_opt:
					print 'found in '+k
				break
		
	else: #fast method
		if organism == 'human':
			r=h5py.File(path_to_use + '/beeswarn_data/h_annot_to_probe.hdf5','r')
		else:
			r=h5py.File(path_to_use + '/beeswarn_data/m_annot_to_probe.hdf5','r')
		
		root= r['root']
		p=None
		for k in root:
			if gene in root[k]:
				p = root[k].get(gene)[...]
			if p != None:
				if not WEB_opt:
					print 'found in '+k
				break
		r.close()
	return p
	pass
	
	
	
#################################################### Functions ############################################


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Plots gene expression accross the hematopoietic system.')
	parser.add_argument('-i','--input', help='Gene 1 to be plotted.')
	parser.add_argument('-j','--input2', help='Gene 2 to be plotted.')#, default='NULL')
	parser.add_argument('-f','--foldchange', action='store_const', const=True, default=False, help='Plot foldchange of gene against normal counterpart.')
	parser.add_argument('data_seq',help='dummy argument - remove when common2 has been updated', type=argparse.FileType('r'), default=sys.stdin )
	parser.add_argument('-o','--organism', type=str,choices= ['human','mouse'] , default='human', help='use human or mouse dataset.'  )
	parser.add_argument('-l','--data', help='Data kind to be plotted.', choices= ['leukemia','normal','both'], default='both')
	parser.add_argument('-lg1','--log2gene1', action='store_const', const=True, default=False, help='Give output of gene 1 in log2 transformed axis')
	parser.add_argument('-lg2','--log2gene2', action='store_const', const=True, default=False, help='Give output of gene 2 in log2 transformed axis')
	args = parser.parse_args()

	in_gene = args.input
#hack for second gene.
	in_gene2 = args.input2
	if in_gene2 == 'NULL':
		in_gene2 = None
	
	log2gene1 = args.log2gene1
	log2gene2 = args.log2gene2

	fold_change = args.foldchange
	organism = args.organism
	data_to_use = args.data
						# CONFIG FILES #
##################################################################################	
##################################################################################	
	#read in configuration file:
#	config = import_sample_data(path_to_use+'JARID2.ini')
	config = import_sample_data('/webdata/servers.binf.ku.dk/programs/shs/JARID2_w.ini') # WEB server.
	
	path_to_use = config[0][1]	#first line is the path
	WEB_opt = int(config[1][1]) #second line is the web option.
##################################################################################	
			
######################################################
#	#sys.path.append('/Users/frederikbagger/BRIC')	 #
######################################################
	
	# load appropriate data file: 
	
	if organism == 'mouse':
		all_data = pickle_load( path_to_use + '/beeswarn_data/nl_mouse_data.pkl')
	else: # for human there are three choices : ['leukemia','normal','both'], if the fc flag is activated, then only leukemia are plotted.
		if not fold_change:
			if data_to_use == 'leukemia' :
				all_data = pickle_load(path_to_use + '/beeswarn_data/all_data_expr.pkl')
			elif data_to_use == 'normal' :
				all_data = pickle_load(path_to_use + '/beeswarn_data/nl_human_data.pkl')
			else:
				merged_data = pickle_load( path_to_use + '/beeswarn_data/nl_human_data.pkl')	   # NOTA : nl_human_data will be updated as soon as the Cancer vs normal 
																								   #		paper is out.
				all_data = pickle_load(path_to_use + '/beeswarn_data/all_data_expr.pkl')
				
				all_data['colnames'] = numpy.concatenate([all_data['colnames'],merged_data['colnames']],axis=0)
				all_data['data'] = numpy.concatenate([all_data['data'],merged_data['data']],axis=1)
				
		else:
			all_data = pickle_load(path_to_use + '/beeswarn_data/all_data_fc.pkl')
#
###########################################################

	# get list of probes for the different genes.
	
	probes_names = find_probes(in_gene, organism)
	
	if probes_names == None:
		print 'Gene name not found on array..'
		sys.exit(2)
	
	probes_to_use=[]
	for g in probes_names:
		try:
			probes_to_use.append(all_data['rownames'].tolist().index(g))
		except Exception,  e:
			pass
	
	if in_gene2 != None:
		
		probes_names2 = find_probes(in_gene2, organism)
		
		if probes_names2 == None:
			print 'Gene 2 name not found on array..'
			sys.exit(2)
			
		probes_to_use2=[]
		for g in probes_names2:
			try:
				probes_to_use2.append(all_data['rownames'].tolist().index(g))
			except Exception,  e:
				pass
		
#########################################################
	
	
	if in_gene2 == None : # just plot the data
#		print probes_to_use
		data=[]
		if len(probes_to_use) == 0:
			print 'Gene not found on array..'
			sys.exit(999)

		a=all_data['data'][probes_to_use,:]
		maxes = numpy.argmax(numpy.abs(a), axis=0)
		for index,i in enumerate(maxes):
			d=all_data['data'][probes_to_use,index][i]
			data.append(d)
			
		
		if organism == 'mouse':
			genename=in_gene.capitalize()
		else: #it is human
			genename= in_gene.upper()
		
		if fold_change:                                    
			beeswarm(numpy.array(data),all_data['colnames'], title = genename, out= in_gene+'_fc.pdf', plot_median=True, use_log=log2gene1) 
			beeswarm(numpy.array(data),all_data['colnames'], title = genename, out= in_gene+'_fc.png', plot_median=True, use_log=log2gene1)
		else:                                              
			beeswarm(numpy.array(data),all_data['colnames'], title = genename, out= in_gene+'.pdf', plot_median=True, use_log=log2gene1)
			beeswarm(numpy.array(data),all_data['colnames'], title = genename, out= in_gene+'.png', plot_median=True, use_log=log2gene1)
                                                           

		
	else:	   # here we plot the correlation
		import matplotlib.pyplot as plt
		fc=[]
		fc2=[]
		all_names=all_data['colnames']
		if not fold_change:
			for i,s in enumerate(all_data['colnames']):
				maxi=0
				for j in all_data['data'][probes_to_use,i]:
					if j>0:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(j))
							signe=numpy.sign(j)
					else:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(-j))
							signe=numpy.sign(j)
				fc.append(maxi*signe)
				maxi=0
				for j in all_data['data'][probes_to_use2,i]:
					if j>0:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(j))
							signe=numpy.sign(j)
					else:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(-j))
							signe=numpy.sign(j)
				fc2.append(maxi*signe)

		else:
			for i,s in enumerate(all_data['colnames']):
				maxi=0
				for j in all_data['data'][probes_to_use,i]:
					if numpy.abs(j) > maxi:
						maxi=numpy.abs(j)
						signe=numpy.sign(j)
				fc.append(maxi*signe)

				maxi=0
				for j in all_data['data'][probes_to_use2,i]:
					if numpy.abs(j) > maxi:
						maxi=numpy.abs(j)
						signe=numpy.sign(j)
				fc2.append(maxi*signe)
#		print 'using %d probeset(s) for %s'%(len(probes_to_use),in_gene)
#		print 'using %d probeset(s) for %s'%(len(probes_to_use2),in_gene2)

		colors=['#ed2921' , '#5f14fa' ,'#fac514', '#0c9ffa' , '#89fa0c' ,  '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914','#ed2921' , '#5f14fa' , '#0c9ffa' , '#89fa0c' , '#fac514' , '#ed1dc5' , '#111dfb' , '#05fb4c' , '#fbea09' , '#ff6914']
		markers = ['s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ]
		
		all_classes_order=[]
		all_classes={}
		for index,classe in enumerate(all_names):
			if classe not in all_classes:
				all_classes_order.append(classe)
				all_classes[classe]=[]
				all_classes[classe].append(index)
			else:
				all_classes[classe].append(index)
	#	print all_classes


		#compute corr coef line.
		lin_fit = numpy.polyfit(numpy.array(fc),numpy.array(fc2),1) # this returns the coef of the polynomial fit
		corr_coef = numpy.corrcoef(numpy.array(fc),numpy.array(fc2))[0][1] # this is R
		line_x = numpy.linspace(numpy.array(fc).min(),numpy.array(fc).max()) # this is to have some points to actually draw the line. 

		plots=[]
		for i,classe in enumerate(all_classes):
			a=plt.plot(numpy.array(fc)[all_classes[classe]],numpy.array(fc2)[all_classes[classe]],'o',alpha=.5, color=colors[i], marker=markers[i], label=classe)
			plots.append(a)
		plots.append( plt.plot(line_x, line_x*lin_fit[0] + lin_fit[1] , '--b', label='$R^2$ = %.2f'%(corr_coef*corr_coef) )) #we append the plot of the line here
		kwarg={'size':6 }
		
		plt.legend(loc='upper right', prop=kwarg)
	
		if log2gene1:
			plt.xscale('log', basex=2)
			plt.xlabel('log2 expression of %s'%in_gene)
			plt.xlim(xmax=plt.xlim()[1]+2^10) #make room for the legend
		else:
			plt.xlabel('%s'%in_gene)
			plt.xlim(xmax=plt.xlim()[1]+1000) #make room for the legend
		
		if log2gene2:
			plt.yscale('log', basey=2)
			plt.ylabel('log2 expression of %s'%in_gene2)
			plt.xlim(xmax=plt.xlim()[1]+2^10) #make room for the legend
		else:
			
			plt.xlabel('%s'%in_gene)
			plt.xlim(xmax=plt.xlim()[1]+1000) #make room for the legend
		
		
		
		if organism == 'mouse':
			genename1=in_gene.capitalize()
			genename2=in_gene2.capitalize()
		else: #it is human
			genename1= in_gene.upper()
			genename2= in_gene2.upper()
		
		if not fold_change:
			plt.title('%s vs %s'%(genename1,genename2))
		else:
			plt.title('%s vs %s (versus normal counterpart)'%(genename1,genename2))


		if not fold_change:
			pylab.savefig('%s.pdf'%in_gene)
			pylab.savefig('%s.png'%in_gene)
		else:
			pylab.savefig('%s_fc.png'%in_gene)
			pylab.savefig('%s_fc.pdf'%in_gene)





##################################################################################	
#MAKE HTML
	head='''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
		<html xmlns="http://www.w3.org/1999/xhtml" lang="en">
		<head>
		<title> Servers.binf.ku.dk | SHS </title>
		<meta name="language" content="en" />
		<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
		<style type="text/css">
		body {
		font-family: 'ArialMT', 'Arial', 'sans-serif';
		font-size: 14px;
		color: #3E3535;
		background-color: #FFFFFF;
		line-height: 20px;
		max-width: 630px;
		margin: 20px auto 20px auto;
		padding: 0px;
		}

		h1.title {font-size: 36px;}
		h2.title {font-size: 36px;}
		pre		 {font-size: 12px;}
		boring	 {font-size: 12px}
		</style>	
		</head>
		<body>
		<br><p style="font-size:150%;font-weight:bold;text-align:center">Welcome to the <span style="color:silver;font-size:250%">HemaExplorer</span></p><p style="text-align:center"><i><FONT COLOR="#000000"> --- See the expression of your favourite gene in the hematopoietic system --- </FONT></i></p>
		<!-- ##### ##### HEADER END ##### ##### -->'''
	
	AMLabr='''<br>Abrieviations:<br>				
		<table border="0">
		<tr>
		<td> AMLI_ETO      </td> <td>  AML with t(8;21)                                  </td>
		</tr>
		<tr>
		<td> APL       </td> <td>  AML with t(15;17)                                    </td>
		</tr>
		<tr>
		<td> AML with inv(16)/t(16;16)  </td> <td>  AML with inv(16)/t(16;16)            </td>
		</tr>
		<tr>
		<td> AML with t(11q23)/MLL   </td> <td>  AML with t(11q23)/MLL                   </td>
		</tr>
		</table>'''
		
	Normabr='''	<br>Abrieviations:<br>				
		<table border="0">
		<tr>
		<td> HSC_BM </td> <td> Hematopoietic stem cells from bone marrow                  </td>
		</tr>
		<tr>
		<td> Early HPC_BM </td> <td> Hematopoietic progenitor cells from bone marrow      </td>
		</tr>
		<tr>
		<td> CMP </td> <td> Common myeloid progenitor cell                                </td>
		</tr>
		<tr>
		<td> GMP</td> <td> Granulocyte macrophage progenitor                             </td>
		</tr>
		<tr>
		<td> MEP </td> <td> Megakaryocyte-erythroid progenitor cell                       </td>
		</tr>
		<tr>
		<td> PM_BM </td> <td> Promyelocyte from bone marrow                               </td>
		</tr>
		<tr>
		<td> MY_BM </td> <td> Myelocyte from bone marrow                                  </td>
		</tr>
		<tr>
		<td> PMN_BM </td> <td> Polymorphonuclear cells from bone marrow                   </td>
		</tr>
		<tr>
		<td> PMN_PB </td> <td> Polymorphonuclear cells from peripheral blood              </td>
		</tr>
		</table>'''
		
	AMLandNormabr='''	<br>Abrieviations:<br>				
		<table border="0">
		<tr>
		<td> HSC_BM </td> <td> Hematopoietic stem cells from bone marrow                  </td>
		</tr>
		<tr>
		<td> Early HPC_BM </td> <td> Hematopoietic progenitor cells from bone marrow      </td>
		</tr>
		<tr>
		<td> CMP </td> <td> Common myeloid progenitor cell                                </td>
		</tr>
		<tr>
		<td> GMP</td> <td> Granulocyte monocyte progenitors                              </td>
		</tr>
		<tr>
		<td> MEP </td> <td> Megakaryocyte-erythroid progenitor cell                       </td>
		</tr>
		<tr>
		<td> PM_BM </td> <td> Promyelocyte from bone marrow                               </td>
		</tr>
		<tr>
		<td> MY_BM </td> <td> myelocyte from bone marrow                                  </td>
		</tr>
		<tr>
		<td> PMN_BM </td> <td> Polymorphonuclear cells from bone marrow                   </td>
		</tr>
		<tr>
		<td> PMN_PB </td> <td> polymorphonuclear cells from peripheral blood              </td>
		</tr>
		<tr>
		<td> AMLI_ETO      </td> <td>  AML with t(8;21)                                  </td>
		</tr>
		<tr>
		<td> APL       </td> <td>  AML with t(15;17)                                    </td>
		</tr>
		<tr>
		<td> AML with inv(16)/t(16;16)  </td> <td>  AML with inv(16)/t(16;16)            </td>
		</tr>
		<tr>
		<td> AML with t(11q23)/MLL   </td> <td>  AML with t(11q23)/MLL                   </td>
		</tr>
		</table>'''

	MouseNormabr='''	<br>Abrieviations:<br>				
		<table border="0">
		<tr>
		<td> LT_HSC </td> <td> Long term Hematopoietic stem cell                  </td>
		</tr>
		<tr>
		<td> ST_HSC </td> <td> Short term Hematopoietic stem cell      </td>
		</tr>
		<tr>
		<td> LMPP </td> <td> Lymphoid-primed multipotential progenitors                   </td>
		</tr>
		<tr>
		<td> CLP</td> <td> Common lymphoid progenitor cells                             </td>
		</tr>
		<tr>
		<td> ETP </td> <td> Early T-cell progenitor                      </td>
		</tr>
		<tr>
		<td> ProB </td> <td> Pro-B cell                               </td>
		</tr>
		<tr>
		<td> PreB </td> <td> Pre-B cell                                  </td>
		</tr>
		<tr>
		<td> IgM+SP </td> <td> Immunoglobulin M positive side population cells                   </td>
		</tr>
		<tr>
		<td> CD4 </td> <td> CD4 cells              </td>
		</tr>
		<tr>
		<td> NKmature </td> <td>  Mature natural killer cells                              </td>
		</tr>
		<tr>
		<td> GMP </td> <td>  Granulocyte monocyte progenitors                            </td>
		</tr>
		<tr>
		<td> MkE </td> <td>  Megakaryocyte erythroid precursors           </td>
		</tr>
		<tr>
		<td> MkP </td> <td> Megakaryocyte precursor                   </td>
		</tr>
		<tr>
		<td> PreCFUE </td> <td> Pre-colony-forming unit erythroid cells                   </td>
		</tr>
		<tr>
		<td> CFUE </td> <td> Colony-forming unit erythroid cells                   </td>
		</tr>
		<tr>
		<td> ProE </td> <td> Erythroid progenitor cells                 </td>
		</tr>
		</table>'''


	singletxt='''	<b>Single gene lookup</b><br>
		Each dot in the plot corresponds the expression of '''+in_gene+''' in a microarray. Horizontal lines represent the median expression for each class of cells. The y-axis is the log2 expression of the gene.'''
	
	singletxt_fc='''	<b>Single gene lookup</b><br>
	Each dot in the plot corresponds the foldchange of '''+in_gene+''' in a microarray. Horizontal lines represent the median expression for each class of cells. The y-axis is the log2 expression of the gene.'''

	if in_gene2 is not None:
		cortxt='''<b>Correlation</b><br>
		Each dot in the plot represent a microarray experiment with different
		cell types, as specified in the legend of the plot. Expression for gene '''+in_gene+''' is given on the x-axis and expression for 
		'''+in_gene2+''' is given on the y-axis. The stippled line represent a theoretical perfect correlation and the R<sup>2</sub> value for the correlation is given in the legend.<br>'''

		cortxt_fc='''<b>Correlation</b><br>
		Each dot in the plot represent a microarray experiment with different
		cell types, as specified in the legend of the plot. The foldchange for '''+in_gene+''' is given on the x-axis and foldchange for '''+in_gene2+''' is given on the y-axis. The stippled line represent a theoretical perfect correlation and the R<sup>2</sub> value for the correlation is given in the legend.<br>'''
		


	
	
	if WEB_opt: # some display
		if in_gene is not None:
			if not fold_change:
				print head
				print '<img src=\"%s.png\" align=center>' % in_gene
				# for human there are three choices : ['leukemia','normal','both'], if the fc flag is activated, then only leukemia are plotted.
				if in_gene2 is None:
					print singletxt
				else:
					print cortxt
				
				if organism == 'mouse':
					print MouseNormabr
				elif data_to_use == 'leukemia' :
					print AMLabr
				elif data_to_use == 'normal' :
					print Normabr
				else:
					print AMLandNormabr 
				print '<br><a target=\"_blank\" href=\"%s.pdf\" title=\"\">Get plot in pdf format</a>' % in_gene
				
			else:
				print head				
				print '<img src=\"%s_fc.png\" align=center>' % in_gene
				if in_gene2 is None:
					print singletxt
				else:
					print cortxt
				if data_to_use == 'leukemia' :
					print AMLabr
				elif data_to_use == 'normal' :
					print Normabr
				else:
					print AMLandNormabr 
				
				
				print '<a target=\"_blank\" href=\"%s_fc.pdf\" title=\"\">Get plot in pdf format</a>' % in_gene
			# print 'gene is #%s#' % in_gene
		else:
			print "Please provide a gene!"
			sys.exit(2)
	else:
		if in_gene is not None:
			print 'gene is #%s#' % in_gene
			if in_gene2 is not None:
				print 'gene 2 is #%s#' % in_gene2
		else:
			print "Please provide a gene!"
