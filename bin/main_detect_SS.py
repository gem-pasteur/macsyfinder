import sys, os

# Script pour la detection des systemes de secretion a partir de profils. 

SYST=['T1SS', 'T2SS', 'T3SS', 'T4SS', 'T5SS', 'T6SS', 'T4P', 'Tad', 'Flagellum']
HMMER="hmmsearch"
RES_SEARCH_DIR="./datatest/res_search"
RES_SEARCH_SUFFIX=".search_hmm.out"
EVALUE_RES=1
PROFILE_DIR="./profiles"

PROFILE_SUFFIX=".fasta-aln_edit.hmm"

RES_EXTRACT_SUFFIX=".res_hmm_extract"

def is_syst_avail(syst):
	
	if(SYST.count(syst)>0):
		return True
	else:
		return False
		
def get_syst_def_file(syst):
	return "DEF/%s.def"%syst

class System:
	'''Classe definissant un SYSTEME a detecter
	
	'''
	
	def __init__(self, system, deffile, profile_dir=PROFILE_DIR, profile_suffix=PROFILE_SUFFIX):
		self.system=system
		if not os.path.exists(deffile):
			raise IOError, "The file %s does not exist"%deffile
		if not os.path.exists(profile_dir):
			raise IOError, "The directory %s does not exist"%profile_dir
			
		self._deffile=deffile
		self._profiledir=profile_dir
		self._profilesuffix=profile_suffix
		(self._listgenes, self._dicogenes)=self.init_genes_from_def_file(self._profiledir, self._profilesuffix)
	
	def __repr__(self):
		format=self.get_format()
		list_genes=self.get_syst_genes()
		
		lines="%s\nDefinition file: %s\nGenes list: %s"%(self.system, self.get_syst_def_file(), self.get_syst_genes())
		
		return lines

	def get_syst_def_file(self):
		return self._deffile
	
	def get_syst_genes(self):
		return self._listgenes

	def get_dico_Genes(self):
		return self._dicogenes
		
	def get_gene_lg(self, gene):
		(self.get_dico_Genes()[gene]).get_profile_lg()
	
	def init_genes_from_def_file(self, profile_dir, profile_suffix):
		genes=[]
		dico_Genes={}
		filename=get_syst_def_file(self.system)
		def_file=open(filename)
		for l in def_file:
			if l!='\n':
				gene=l[:-1].strip()
				if not dico_Genes.has_key(gene):
					genes.append(gene)
					profile="%s/%s%s"%(profile_dir, gene, profile_suffix)
					dico_Genes[gene]=Gene(gene, profile, self.system)
					
		def_file.close()
	
		return (genes, dico_Genes)
		
	def detect_hmmer_genes(self, hmmer_exe, sequence_db, evalue_thresh, outfile_dir, outfile_suffix, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix, data_type=""):
		allGenes=self.get_dico_Genes()
		outfiles=[]
		for el in self.get_syst_genes():
			gene=allGenes[el]
			
				
			#print "\nPerforming Hmmer search for :"
			print(gene)
			
			search=HmmerGeneSearchPerformer(hmmer_exe, gene, sequence_db, evalue_thresh, outfile_dir, outfile_suffix)
			search.perform_search()
			
			#print "- outfile: %s"%search._outfile
			print "- outfile: %s"%search.get_outfile()
			
			#res=ResGeneSearchExtractor(search, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix)
			if (data_type=="gembase"):
				res=ResGeneSearchExtractorGembase(search, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix)
			
			else:
				# NEW	
				res=ResGeneSearchExtractorUnorderedDataset(search, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix)
				#res.extract_res()
			print "Hits extracted to: %s"%res.outfile
			outfiles.append(res.outfile)
		return outfiles
			#break
			
		


class Gene:
	'''Classe definissant un GENE a detecter
	
	'''
	
	def __init__(self, name, profile, system):
		if not os.path.exists(profile):
			raise IOError, "The file %s does not exist"%profile
			
		self.name=name
		self.system=system
		self._profile=profile
		self._profile_lg=self.init_lg_profile(self._profile)
	
	def __repr__(self):
		lines="---------\n- %s -\n---------\n- profile: %s\n- length: %d"%(self.get_gene_name(), self.get_profile(), self.get_profile_lg())
		
		return lines

	def init_lg_profile(self, profile):
		fich=open(profile)
		for l in fich:
			if l.startswith("LENG"):
				length=int(l.split()[1])
				break
		fich.close()		
		return length

	def get_gene_name(self):	
		return self.name
		
	def get_profile(self):
		return self._profile

	def get_profile_lg(self):
		return self._profile_lg
	
	
class HmmerGeneSearchPerformer:
	'''Classe definissant une recherche Hmmer pour: 
		- un GENE donne et 
		- une BD de sequences	
	'''
	
	def __init__(self, hmmer_exe, gene_obj, sequence_db, evalue_thresh, outfile_dir, outfile_suffix):
		self.gene=gene_obj
		self._geneprofile=self.gene.get_profile()
		self._sequence_db=sequence_db
		#self._profile_lg=self.init_lg_profile(self._profile)
		self._outfile="%s/%s%s"%(outfile_dir, self.gene.get_gene_name(), outfile_suffix)
		self._evalue_thresh=evalue_thresh
		self._hmmer_exe=hmmer_exe
		
		#self.perform_search(self._hmmer_exe)
	
	def perform_search(self):
	
		if not os.path.exists(self._outfile):
			res=os.system("%s -o %s -E %d %s %s"%(self._hmmer_exe, self._outfile, self._evalue_thresh, self._geneprofile, self._sequence_db))
		else:
			res=0 # TMp
		
		if (res==0):
			print "Hmmer search done. "
			#return outfile
		else:
			print "Hmmer search failed."
			#return ""		
	
	def get_outfile(self):
		return self._outfile
		
	def get_evalue(self):
		return self._evalue_thresh
		
	def get_hmmer_exe(self):
		return self._hmmer_exe
	
class ResGeneSearchExtractor:
	'''Classe definissant une extraction de resultats d'une recherche Hmmer pour: 
		- une HmmerGeneSearch donnee
		- un seuil de i-evalue donne
		- un seuil de %coverage donne		
	'''
	
	def __init__(self, hmmergenesearch, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix):
		self._i_evalue=i_evalue_thresh
		self._pc_cov=pc_length_coverage_thresh
		self.outfile="%s/%s%s"%(resextract_dir, hmmergenesearch.gene.get_gene_name(), resextract_suffix)
		evalue=hmmergenesearch.get_evalue()
		self.profile_lg=hmmergenesearch.gene.get_profile_lg()
		
		if(self._i_evalue> evalue):
			print "i-evalue overriden by e-value threshold of HMMER search: %f"%evalue
			self._i_evalue=evalue
		


class ResGeneSearchExtractorUnorderedDataset(ResGeneSearchExtractor):
	'''Classe definissant une extraction de resultats d'une recherche Hmmer pour: 
		- une HmmerGeneSearch donnee
		- un seuil de i-evalue donne
		- un seuil de %coverage donne
		- des donnees non ordonnees !		
	'''
	
	def __init__(self, hmmergenesearch, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix):
		# Appel au constructeur de la super-classe
		ResGeneSearchExtractor.__init__(self, hmmergenesearch, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix)
		#self.extract_res(hmmergenesearch.get_outfile(),self.outfile, hmmergenesearch.gene.get_gene_name(), hmmergenesearch.gene.system, i_evalue_thresh, profile_lg, pc_length_coverage_thresh)
		self.extract_res(hmmergenesearch.get_outfile(),self.outfile, hmmergenesearch.gene.get_gene_name(), hmmergenesearch.gene.system, i_evalue_thresh, self.profile_lg, pc_length_coverage_thresh)


	def extract_res(self, search_outfile, res_outfile, gene_name, gene_system, i_evalue_thresh, gene_profile_lg, pc_length_coverage_thresh):
		fich=open(search_outfile)
		lines=fich.readlines()
		fich.close()
		outlines=[]
		i=0
		while i<len(lines)-4:
		
			line=lines[i]
			if(not line.startswith(">> ")):
				i+=1
			else:
				hit_id=line.split()[1]
				#print hit_id
				#line=fich.readline()	
				#line=fich.readline()
				if(lines[i+2].startswith(" ---   ------")):
					i+=3
					line=lines[i]
					#line=lines
					cur_lg=0
					cur_i_eval=""
					while(len(line)>2 and line.split()[0]!="Alignments"):
						fields=line.split()
						#print fields
						#print i_evalue_thresh
						#print "%f %f"%(float(fields[4]), i_evalue_thresh)
						if(len(fields)>1 and float(fields[5]) <= i_evalue_thresh):
							cov=(float(fields[7])-float(fields[6])+1)/gene_profile_lg
							if (cov >= pc_length_coverage_thresh):
								i_eval=fields[5]
								score=fields[2]
								#print hit_id
								#print "%s\t%f\t%f\t%f"%(hit_id, float(fields[4]), i_evalue_thresh, cov)
								outlines.append("%s\t%s\t%s\t%s\t%s\t%f\n"%(hit_id, gene_name, gene_system, i_eval, score, cov))
								#print cov
						i+=1
						line=lines[i]	
				else:
					i+=1
				#break
			#line=fich.readline()	
			#i+=1
			
		fich=open(res_outfile, "w")
		fich.writelines(outlines)
		fich.close()


class ResGeneSearchExtractorGembase(ResGeneSearchExtractor):
	'''Classe definissant une extraction de resultats d'une recherche Hmmer specifique a Gembase (nom replicon et position inclus dans le nom du gene) pour: 
		- une HmmerGeneSearch donnee
		- un seuil de i-evalue donne
		- un seuil de %coverage donne		
	'''
	

	def __init__(self, hmmergenesearch, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix):
		# Appel au constructeur de la super-classe
		ResGeneSearchExtractor.__init__(self, hmmergenesearch, i_evalue_thresh, pc_length_coverage_thresh, resextract_dir, resextract_suffix)
		#self.extract_res(hmmergenesearch.get_outfile(),self.outfile, hmmergenesearch.gene.get_gene_name(), hmmergenesearch.gene.system, i_evalue_thresh, profile_lg, pc_length_coverage_thresh)
		self.extract_res(hmmergenesearch.get_outfile(),self.outfile, hmmergenesearch.gene.get_gene_name(), hmmergenesearch.gene.system, i_evalue_thresh, self.profile_lg, pc_length_coverage_thresh)
		
	def extract_res(self, search_outfile, res_outfile, gene_name, gene_system, i_evalue_thresh, gene_profile_lg, pc_length_coverage_thresh):
		fich=open(search_outfile)
		lines=fich.readlines()
		fich.close()
		outlines=[]
		i=0
		while i<len(lines)-4:
		
			line=lines[i]
			if(not line.startswith(">> ")):
				i+=1
			else:
				hit_id=line.split()[1]
				#print hit_id
				#line=fich.readline()	
				#line=fich.readline()
				if(lines[i+2].startswith(" ---   ------")):
					i+=3
					line=lines[i]
					#line=lines
					cur_lg=0
					cur_i_eval=""
					while(len(line)>2 and line.split()[0]!="Alignments"):
						fields=line.split()
						#print fields
						#print i_evalue_thresh
						#print "%f %f"%(float(fields[4]), i_evalue_thresh)
						if(len(fields)>1 and float(fields[5]) <= i_evalue_thresh):
							cov=(float(fields[7])-float(fields[6])+1)/gene_profile_lg
							if (cov >= pc_length_coverage_thresh):
								i_eval=fields[5]
								score=fields[2]
								#print hit_id
								#print "%s\t%f\t%f\t%f"%(hit_id, float(fields[4]), i_evalue_thresh, cov)
								
								fields_hit=hit_id.split('_')
								replicon_name=fields_hit[0]
								position_hit=int(fields_hit[1])/10
								#outlines.append("%s\t%s\t%s\t%s\t%s\t%f\n"%(hit_id, gene_name, gene_system, i_eval, score, cov))
								outlines.append("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%f\n"%(hit_id, replicon_name, position_hit, gene_name, gene_system, i_eval, score, cov))
								#print cov
						i+=1
						line=lines[i]	
				else:
					i+=1
				#break
			#line=fich.readline()	
			#i+=1
			
		fich=open(res_outfile, "w")
		fich.writelines(outlines)
		fich.close()

def select_best_hit_same_id(infile, outfile, field_id, field_comp, field_gene, criterion_min_eval=True):
	"""
		Select the best hit when multiple hits for a same ID field. 
	"""
	
	fich=open(infile)
	ficout=open(outfile,"w")
	identif_prev=""
	evalue_prev=float(-1)
	ligne_prev=""
	for l in fich:
		fields=l.split()
	
		#if len(fields)==6:
		identif=fields[field_id]
		#print identif
		evalue=float(fields[field_comp])
		#print evalue
		if identif==identif_prev:
			
			if criterion_min_eval:
				if evalue < evalue_prev:
					ligne_prev=l
					evalue_prev=evalue
			else:
				if evalue > evalue_prev:
					ligne_prev=l
					evalue_prev=evalue
				
		else:
			if ligne_prev!="": # A remettre
				ficout.write(ligne_prev)
			identif_prev=identif
			evalue_prev=evalue
			ligne_prev=l
		

	ficout.write(ligne_prev)
	fich.close()
	ficout.close()


def study_cluster_ordered_data(infile, outfile, max_spacer, col_replicon_name, col_gene_position, col_hit):
	"""
		Cover the list of hits and determine 'patches' of hits, i.e sets of hits non separated by more than 'max_spacer' genes. 
		The input must have: 
			- a col with the replicon name (# of col specified by 'col_replicon_name'), 
			- a col with the gene position (# of col specified by 'col_gene_position'), 
			- a column with the hit (# of col specified by 'col_hit').
	"""
	print "To be implemented !! and to replace study_cluster_gembase???"
	

#def study_cluster_gembase(infile, outfile, diff_pos, col_gene_name=1, col_hit=2):
def study_cluster_gembase(infile, outfile, diff_pos, col_gene_name=1, col_hit=4, col_replicon_name=2, col_gene_position=3):
	"""
		#A remplacer par study_cluster_ordered_data (les donnees de gembase sont ordonnees!)
		
		Cover the list of hits and determine 'patches' of hits, i.e sets of hits non separated by more than 'max_spacer' genes. 
		The input must have: 
			- a col with the gene name (# of col specified by 'col_replicon_name'), 			
			- a column with the gene matched with Hmmer (# of col specified by 'col_hit').
			- a column with the replicon name? 
			- a column with the gene position along the replicon? 
			
	"""
	resfich=open(infile) # Par exemple "all_res.sort_chr"
	lines=resfich.readlines()
	resfich.close()
		
	chr_list=[]
	prev_gene_list=[] # paquet 
	prev_geneid_list=[] # paquet # new
	prev_gene_pos_list=[] # paquet


	for l in lines:
		fields=l.split()
		cur_chr=fields[col_replicon_name-1]
		cur_pos=int(fields[col_gene_position-1])
		cur_gene=fields[col_hit-1]
		cur_geneid=fields[col_gene_name-1] # new
		
		if chr_list.count(cur_chr)==0:
		
			# Traitement du chromosome precedent : 
			if len(chr_list)>0:
				# On peut plus agglomerer, on traite la liste en cours du chr precedent
				if len(chr_gene_list)>=2:
					# On garde la liste en cours ssi elle contient au moins 2 genes
					cur_list_gene_list.append(chr_gene_list)
					cur_list_pos_list.append(chr_gene_pos_list)				
					cur_list_geneid_list.append(chr_geneid_list) # New
					
				# On fait une svg pour le chr precedent... 
				# NEW clusters "relaches" : Ajout d'une liste de listes pour clusters "relaches" et d'un compteur global des core genes dans des clusters "relaches"
				prev_list_gene_list=cur_list_gene_list
				prev_list_geneid_list=cur_list_geneid_list # New
				prev_list_pos_list=cur_list_pos_list
				prev_list_core_genes_nr=cur_list_core_genes_nr	

				##### Traitement du precedent avant de reinitialiser...#####
				nb_core=len(prev_list_core_genes_nr)
				if nb_core>0:
					# On a pour le chr precedent tous les core genes en paquets d'au moins 2 genes:
				
					print "\n***** %s *****"%prev_chr	
					for (liste_gene, liste_pos, liste_geneid) in zip(prev_list_gene_list, prev_list_pos_list, prev_list_geneid_list): # New
						
						print "- paquet -"
						print liste_geneid # New
						print liste_gene
						print liste_pos
						
			# Traitement du nouveau chr : 
			chr_list.append(cur_chr)
			chr_gene_list=[] # paquet 
			chr_gene_list.append(cur_gene)
			
			# NEW:
			chr_geneid_list=[] # paquet 
			chr_geneid_list.append(cur_geneid)
			
			chr_gene_pos_list=[] # paquet 
			chr_gene_pos_list.append(cur_pos)
			
			# Ajout d'une liste de listes pour clusters "relaches" et d'un compteur global des core genes dans des clusters "relaches"
			cur_list_gene_list=[] # liste de paquets # genes matches
			cur_list_geneid_list=[] # liste de paquets # nom gene # NEw
			cur_list_pos_list=[] # liste de paquets
			cur_list_core_genes_nr=[] # liste non redondante de "core genes"
		
		
		else:
				
			# Chromosome deja en cours de traitement... 
			if(cur_pos-prev_pos) <= diff_pos:
				# On agglomere 
				chr_gene_list.append(cur_gene)
				chr_geneid_list.append(cur_geneid) # new
				chr_gene_pos_list.append(cur_pos)
				
				# La liste fait forcement au moins 2 elements... 
				if cur_list_core_genes_nr.count(cur_gene)==0:
					cur_list_core_genes_nr.append(cur_gene)
				if len(chr_gene_list)==2:
					if cur_list_core_genes_nr.count(prev_gene)==0:
						cur_list_core_genes_nr.append(prev_gene)
			else:
				# On peut pas agglomerer
				if len(chr_gene_list)>=2:
					# On garde la liste en cours ssi elle contient au moins 2 genes
					cur_list_gene_list.append(chr_gene_list)
					cur_list_geneid_list.append(chr_geneid_list) # New
					cur_list_pos_list.append(chr_gene_pos_list)

				chr_gene_list=[]
				chr_gene_list.append(cur_gene)	
			
				# NEW:
				chr_geneid_list=[] # paquet 
				chr_geneid_list.append(cur_geneid)
			
				chr_gene_pos_list=[]
				chr_gene_pos_list.append(cur_pos)


			
		prev_chr=cur_chr
		prev_pos=cur_pos
		prev_gene=cur_gene
		prev_geneid=cur_geneid # New
		
		prev_list_core_genes_nr=cur_list_core_genes_nr
		prev_list_gene_list=cur_list_gene_list
		prev_list_geneid_list=cur_list_geneid_list # New
		prev_list_pos_list=cur_list_pos_list
		
		
	###### Traitement du dernier chr !! ######
	nb_core=len(prev_list_core_genes_nr)
	if nb_core>0:
		print "***** %s *****"%prev_chr
		#for (liste_gene, liste_pos) in zip(prev_list_gene_list, prev_list_pos_list):
		for (liste_gene, liste_pos, liste_geneid) in zip(prev_list_gene_list, prev_list_pos_list, prev_list_geneid_list): # nEW
	
			print liste_geneid # NEw
			print liste_gene
			print liste_pos
		
		print "\n"
	

		
def main():

        if len(sys.argv)<3:
                print "Usage: python main_detect_SS.py sequence_db SYST1 SYST2 ..."
                return 1
        else:

		# Recuperation de la base de sequences a interroger
		sequence_file = sys.argv[1]
		
		# Recuperation de la liste des systemes
		systems_list=[]
		dico_syst={}
		for el in sys.argv[2:]:
			if(is_syst_avail(el)):
				systems_list.append(el)
				if not dico_syst.has_key(el):
					dico_syst[el]=System(el, get_syst_def_file(el), PROFILE_DIR, PROFILE_SUFFIX)
				print "%s available"%el
			else:
				print "%s not available"%el
			
		outfiles=[]
		for syst in systems_list:
			print "\nPerforming HMMER search for %s\n"%syst
			#dico_syst[syst].search_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX)
			#outies=dico_syst[syst].detect_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX, 0.001, 0.5, RES_SEARCH_DIR, RES_EXTRACT_SUFFIX)
			#outies=dico_syst[syst].detect_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX, 0.001, 0.3, RES_SEARCH_DIR, RES_EXTRACT_SUFFIX)
			#outies=dico_syst[syst].detect_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX, 0.001, 0.3, RES_SEARCH_DIR, RES_EXTRACT_SUFFIX, "gembase")
			#outies=dico_syst[syst].detect_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX, 0.001, 0.5, RES_SEARCH_DIR, RES_EXTRACT_SUFFIX, "gembase")
			
			# Donnees Gembase
			outies=dico_syst[syst].detect_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX, 0.00001, 0.5, RES_SEARCH_DIR, RES_EXTRACT_SUFFIX, "gembase")
			# TEST metagenomes
			#outies=dico_syst[syst].detect_hmmer_genes(HMMER, sequence_file, EVALUE_RES, RES_SEARCH_DIR, RES_SEARCH_SUFFIX, 0.00001, 0.5, RES_SEARCH_DIR, RES_EXTRACT_SUFFIX)
			
			#print out
			for out in outies:
				outfiles.append(out)
			
		# Recuperation des hits et creation d'un fichier global de resultats? 
		# Besoin de recuperer la taille du profil... Dans 
		print "\nBuilding a joint resfile"
		alloutlines=[]
		for out in outfiles:
			print out
			fich=open(out)
			lines=fich.readlines()
			fich.close()
			for l in lines:
				alloutlines.append(l)
		
		#print alloutlines	
		alloutlines.sort()
		#print alloutlines
		
		# TMP
		infile="datatest/joint_res_file"
		#infile="joint_res_file_Metag"
		
		fich=open(infile, "w")
		fich.writelines(alloutlines)
		fich.close()
		
		outfile="datatest/joint_res_file.sel"
		#outfile="joint_res_file.sel_Metag"
		
		select_best_hit_same_id(infile, outfile, 0, 5, 1)
		
		max_gene_spacer=15 # All except T1SS/T5SS/Flagellum
		#max_gene_spacer=5 # For T1SS
		#max_gene_spacer=25 # For flagellum 
		outclusters="res_file_clusters_%d"%max_gene_spacer
		#study_cluster_gembase(outfile, outclusters, max_gene_spacer, 1, 2)
		
		# A remettre pour etude dans GEmbase
		study_cluster_gembase(outfile, outclusters, max_gene_spacer)
		
		# A implementer : plus general !! 
		#study_cluster_ordered_data(infile, outfile, max_spacer, col_replicon_name, col_gene_position, col_hit)
		
		
		
			
main()	
	
	
	
