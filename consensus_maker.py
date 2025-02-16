from sys import getsizeof
import random
import os
import argparse
import time

var_dct={}
chr_var_coords={}
sequences={}
cons_dct={}
chrom=''
chromosome_list=[]
chrom_starts={}
cons_starts_dct={}
sample_start=None
samples=[]

start_time = time.time()

def chr_format_unify(chr):
    if chr.startswith('chr') == False:
        return 'chr'+chr
    else: 
        return chr
    
def get_elapsed_time():
    return round(time.time()-start_time,3)

parser = argparse.ArgumentParser(
                prog='consensus_maker.py',
                description='''This tool makes consensuses from an input FASTA file and a VCF file containing variant calls.''',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('input_fasta',
                    help='Input FASTA file path.')
parser.add_argument('vcf',
                    help='Input VCF file path.')
parser.add_argument('output_fasta',
                    help='Output FASTA path.')

parser.add_argument('-c','--cons_count', type=int,
                    help='Number of consensuses', default=1000)
parser.add_argument('-l','--cons_length', type=int,
                    help='Length of consensuses',default=160)
parser.add_argument('-m','--min_freq', type=float,
                    help='Minimum allele frequency to include a variant in the consensus',default=0.0)
parser.add_argument('-s','--sorted', action='store_true',help='Sort consensuses by coordinate')
parser.add_argument('-e','--cons_extra_bases', type=int,
                    help='(Advanced) Extra bases (max) to be generated in each consensus to compensate for deletions',default=200)
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()
cons_count=args.cons_count
cons_length=args.cons_length
fasta=args.input_fasta
vcf=args.vcf
output_fasta=args.output_fasta
min_freq=args.min_freq
cons_extra_bases=args.cons_extra_bases
sorted=args.sorted

print('Generating',cons_count,'consensuses of length',cons_length)
print('Minimum allele frequency:',min_freq)

##Parsing VCF file...
print(f'[{get_elapsed_time()}s] Processing VCF file...')
with open(vcf, 'r') as vcf_file:
    for line in vcf_file:
        
        ##Skipping useless headers...
        if line.startswith('##'):
            continue

        ##Recording useful headers...
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            continue

        ##Parsing data row...
        if headers:
            fields=line.strip().split('\t')
            fields[0]=chr_format_unify(fields[0])

            ##Adding useful fields...
            sample_genotype_list=()
            for i in range(len(headers)-9): ##Iterating the indices of samples...
                sample_genotype_list=(*sample_genotype_list,fields[i+9])

            ##Adding the current VCF line to the variant dictionary - (Chr, coord):(REF,ALT,FORMAT,Sample1,Sample2...)    
            var_dct[(fields[0],int(fields[1]))]=(fields[3],fields[4],fields[7],fields[8],*sample_genotype_list)

            ##Adding variant coordinate to a dictionary of variant coordinates - Chr : [coord1, coord2, coord3...]
            if fields[0] not in chr_var_coords.keys():
                chr_var_coords[fields[0]]=[int(fields[1])]
            else:
                chr_var_coords[fields[0]].append(int(fields[1]))
                           
##Finding sample information in VCF headers...
for element in headers:
    if element == 'FORMAT':
        sample_start=True
        continue
    if sample_start:
        samples.append(element)
print(f'[{get_elapsed_time()}s] Samples:',samples)

##Parsing FASTA file...
print(f'[{get_elapsed_time()}s] Processing FASTA file (can take a while)...')
with open(fasta, 'r') as reference:
    for line in reference:
        if line.startswith('>'):
            if chrom:
                sequences[chrom]=seq  
        
            ##Parsing chromosome...
            header_fields=line[1:].strip().split(' ')           
            if 'chromosome' in  header_fields:
                chrom = chr_format_unify(header_fields[header_fields.index('chromosome')+1].strip(','))
            else: 
                chrom=header_fields[0]           
            seq=''

            ##Parsing start position if header contains  (e.g. NCBI format)... 
            if 'ref' in header_fields[0]:
                fasta_range=header_fields[0].split(':')[1]
                fasta_start=fasta_range[0:fasta_range.find('-')]
            else: fasta_start = 0
            chrom_starts[chrom]=int(fasta_start)

        else:
            seq+=line.strip()
    sequences[chrom]=seq ##Parsing last line

##Making consensuses...
print(f"[{get_elapsed_time()}s] Making consensuses...")
while len(list(cons_dct)) < cons_count:

    ##Randomly defining consensus parameters...
    chrom=random.choice(list(sequences.keys())) ##Choosing a random chromosome...
    sample_ID=random.randint(1,len(samples)) ##Choosing a random sample...
    cons_start=random.randint(0,len(sequences[chrom])-cons_length-cons_extra_bases) ##Choosing a consensus start coord...

    #Adding consensus coordinates to a dictionary...
    if chrom in cons_starts_dct.keys():
        cons_starts_dct[chrom].append(cons_start+chrom_starts[chrom])  
    else: 
        cons_starts_dct[chrom] = [cons_start+chrom_starts[chrom]]

    ##Extra bases are needed to compensate for deletions. They get removed later.
    consensus=sequences[chrom][cons_start:cons_start+cons_length+cons_extra_bases]

    ##Making a list of bases in the consensus...
    string_list=list(consensus)

    ##Setting some placeholder values...
    variant_counter=0
    ploidy=0
    cons_chromosome_ID=1
    new_string_list=string_list.copy()
    symbolic=False

    ##Modifying the consensus with variants from the VCF file...
    if chrom in list(chr_var_coords.keys()): ##Checking if the consensus original chromosome is in VCF file...

        ##Checking ploidy of a selected chromosome...
        coord=random.choice(chr_var_coords[chrom])
        ploidy=len(var_dct[chrom, coord][3+sample_ID].split(':')[0].replace('/','|').strip('|').split('|'))

        #Randomly choosing a chromosome for the consensus...
        cons_chromosome_ID=random.randint(1,ploidy)
        for coord in chr_var_coords[chrom]:
            if cons_start<=coord-chrom_starts[chrom]<cons_start+cons_length+cons_extra_bases:
                
                sample_genotype=var_dct[chrom, coord][3+sample_ID].split(':')[0].replace('/','|').strip('|').split('|')[cons_chromosome_ID-1] ##Determining consensus genotype...
                ref_len=len(var_dct[chrom, coord][0]) ##Length of a reference allele in the VCF file

                ##Extracting allele frequency from var_dct...
                info_field=var_dct[chrom,coord][2].split(';')
                allele_frequency=float(next(i for i in info_field if 'AF=' in i)[3:])
              
                ##Inserting ALT allele if conditions are right...
                if sample_genotype != '0' and sample_genotype != '.' and allele_frequency >= min_freq:

                    ##Incrementing the variant counter (can make mistakes if long deletions are present)...
                    if coord-chrom_starts[chrom]<cons_start+cons_length:
                        variant_counter+=1

                    ##Choosing ALT variant...
                    if ',' not in var_dct[chrom, coord][1]:
                        alt_variant = var_dct[chrom, coord][1]
                    else:
                        alt_variant= var_dct[chrom, coord][1].split(',')[int(sample_genotype)-1]

                    if '<' in alt_variant: ##Checking if alt is a structural variant...
                        symbolic= True
                        svlen=int(next(i for i in info_field if 'SVLEN=' in i)[6:])
                    else:        
                        new_string_list[coord-cons_start-chrom_starts[chrom]]=alt_variant ##Inserting simple variant...    

                    ##Deleting simple deletions...
                    if ref_len > 1:                         
                        for i in range(ref_len-1):
                            new_string_list[min(coord-chrom_starts[chrom]-cons_start+i,cons_length+cons_extra_bases-1)]=''

                    ##Deleting symbolic deletions... 
                    if 'DEL' in alt_variant:
                        for i in range(svlen):
                            new_string_list[min(coord-chrom_starts[chrom]-cons_start+i,cons_length+cons_extra_bases-1)]=''
              
    if new_string_list!=string_list:
        consensus=''.join(new_string_list)
        
    cons_dct[(chrom,cons_start)]=(consensus[0:cons_length], variant_counter, sample_ID, ploidy, cons_chromosome_ID)

##Making a list of chromosomes in input FASTA...
for key in cons_starts_dct.keys():
    chromosome_list.append(key)
chromosome_list.sort()   

##Sorting chromosomes and consensus coordinates lists...
if sorted == True:
    print(f'[{get_elapsed_time()}s] Sorting consensuses...')
    for value in cons_starts_dct.values():
        value.sort()
else: print(f'[{get_elapsed_time()}s] Not sorting consensuses.')

print(f'[{get_elapsed_time()}s] Consensuses made, writing them to output...')

##Writing consensuses to output FASTA...
with open(output_fasta, 'w') as output_file:
    for chrom in chromosome_list:
        for coord in cons_starts_dct[chrom]:
            cons=cons_dct[(chrom, coord-chrom_starts[chrom])]

            ##Writing header...
            output_file.write(f">{chrom}, sample '{samples[cons[2]-1]}', {coord}-{coord+cons_length} bp, ploidy = {cons[3]}, chromosome ID #{cons[4]}, {cons[1]} alt. variants \n")
            
            ##Writing sequence...
            fasta_line_length=0
            while fasta_line_length < cons_length:
                output_file.write(cons[0][fasta_line_length:min(fasta_line_length+80,len(cons[0]))]+'\n')
                fasta_line_length+=80

fasta_size=os.stat('consensuses.fasta').st_size
print(f"[{get_elapsed_time()}s] Consensus FASTA created, size: {len(cons_dct)} sequences, {int(fasta_size/1000)} KB.")