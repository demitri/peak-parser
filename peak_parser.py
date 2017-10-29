#! /usr/bin/env python3

#object will be path defined by user when this method is called:
#path = sys.argv[1]
#newfile = Bedfile(path=path)
class Bedfile(object):
	'''
	Class to read a bedfile type file.
	'''
	def __init__(self, path=None): #standard idiom when defining classes

		self.chromosomes = dict() # key=chromosome name, value=Chromosome object
		self.filepath = path
				
		with open(path) as bedfile:
			for line in bedfile:
				line = line.rstrip()
				
				fields = line.split("\t")
				chr_name = fields[0]

				# get existing or create new chromosome				
				if chr_name in self.chromosomes:
					chromosome =  self.chromosomes[chr_name] #chromosome_name is associated with the chromosome value, which in the else statement below we turn into a Chromosome object list
				else:
					chromosome = Chromosome(name=chr_name)
					self.chromosomes[chr_name] = chromosome

				# create new peak from row
				new_peak = Peak() #create new object with Peak class method
				new_peak.start = int(fields[1]) #convert the string to an int
				new_peak.end = int(fields[2])
				new_peak.signal = float(fields[6])
				peak_average = str(fields[3]).split(":")
				new_peak.mean = peak_average[1]
				#append the list of gene properties to the chromosome dictionary
				chromosome.peaks.append(new_peak) #chromosome is an object containing a list of objects
			#print(chromosome.peaks)
	
	def chromosome_names(self):		
		return self.chromosomes.keys()  	
							
class Chromosome(object): #this is an ontology class designation object: we want to be able to attribute object lists to our dictionaries. Making a list a type Chromosome class allows us include properties in a list form.
	def __init__(self,  name=None):
		self.name = name
		self.length = None	# number of base pairs
		self.peaks = list() # will be used for a list of peak objects
		self.genes = list() # will be used for a list of gene objects
		self.degenes = list() # will be used for a list of differentially expressed gene objects---RNAseq foldchange file
	
class Peak(object):
	def __init__(self):

		#parenthesis .start() means I'm calling a function. .start means it is a property (variable)
		#properties for holding peak information. Must be defined here as something or nothing, in this case, None.
		#avoid placeholder properties until you know what you want: it takes up memory and other things
		self.start = None
		self.end = None
		self.signal = None
		self.mean = None
		
	def width(self):
		return (self.end - self.start) #must use this variable (property) within this class because it hasn't been defined elsewhere.
	
	def intersects_with(self, other_peak):
		pass # return True or False

		
class Gff_file(object):
	def __init__(self, path=None):
		
		self.chromosomes= dict()
		self.filepath = path
		#parse file as normal				
		with open(path) as gf:
			for line in gf:
				line = line.rstrip()
				fields = line.split("\t")
				chr_name = fields[0] #setting chromosome name
			# get existing or create new chromosome				
				if chr_name in self.chromosomes:
					chromosome =  self.chromosomes[chr_name] #add chromosome name as a key, this will NOT replace the values already associated with an existing identical key. Else will turn chromosome into a Chromosome object
				else:
					chromosome = Chromosome(name=chr_name) #creates Chromosome object with key name as the chromosome name
					self.chromosomes[chr_name] = chromosome #fills the dictionary with the Chromosome object.
				#get line only if it is a gene
				if "gene" == fields[2]:
					gene = Gene() #set new class so we can attribute properties.
					if "+" == fields[6]: #identify if TSS is for a gene on positive or negative strand
						tss = int(fields[3])
						gene.tss = tss
						gene.tss_1kb_upstream = (gene.tss - 1000)
					elif "-" == fields[6]:
						tss = int(fields[4]) # 5' to 3' direction changes on opposite strand
						gene.tss = tss
						gene.tss_1kb_upstream = (gene.tss + 1000) 				
					gene.start = int(fields[3]) #start is tss for positive strand -1000 for peak analysis
					gene.stop = int(fields[4]) #stop is tss for negative strange +1000 for peak analysis
					gene.strand = fields[6]
					geneid = fields[8].split(";")
					gene.ID = geneid[0]
					#append the list of gene properties to the chromosome dictionary
					chromosome.genes.append(gene)
				
class Gene(object):
	def __init__(self):
	
		self.start = None
		self.stop = None
		self.ID = None
		self.strand = None
		self.tss = None
		self.tss_1kb_upstream = None

class DE_genes(object):
	def __init__(self, path=None):
		self.expression_values_dict = dict()
	
		with open(path) as defile:
			for line in defile:
				line = line.rstrip()
				fields = line.split('\t')
				gene_name = fields[2]

				if gene_name not in self.expression_values_dict:
					gene_value_list_object = Chromosome(name = gene_name) #calling Chromosome class here to populate list of objects later					
					self.expression_values_dict[gene_name] = gene_value_list_object
					gene = DEexpression()
					gene.cluster = fields[0]
					gene.ID = fields[2]
					gene.foldchange = fields[3:len(fields)-1]
					gene.annotation = fields[-1]
					
				gene_value_list_object.degenes.append(gene)
		
			
		
class DEexpression(object):
	def __init__(self):
	
		self.cluster = None
		self.ID = None
		self.foldchange = None
		self.annotation = None
	
	
	
	
	