#Tom's Script for Editing Manta Files Into Something DupHold Can use

import sys
import cyvcf2
from cyvcf2 import Writer
import numpy as np

vcf = cyvcf2.VCF(sys.argv[1])
output_vcf = sys.argv[2]

vcf.add_info_to_header({'ID': 'STRANDS', 'Description': 'Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+)', 'Type': 'String', 'Number': '.'})
vcf.add_info_to_header({'ID': 'PRPOS', 'Description': 'Breakpoint probability dist', 'Type': 'String', 'Number': '1'})
vcf.add_info_to_header({'ID': 'PREND', 'Description': 'Breakpoint probability dist', 'Type': 'String', 'Number': '1'})
vcf.add_info_to_header({'ID': 'STRANDS', 'Description': 'Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--', 'Type': 'String', 'Number': '.'})
vcf.add_info_to_header({'ID': 'SU', 'Description': 'Number of pieces of evidence supporting the variant across all samples', 'Type': 'Integer', 'Number': '.'})
vcf.add_info_to_header({'ID': 'PE', 'Description': 'Number of paired-end reads supporting the variant across all samples', 'Type': 'Integer', 'Number': '.'})
vcf.add_info_to_header({'ID': 'SR', 'Description': 'Number of split reads supporting the variant across all samples', 'Type': 'Integer', 'Number': '.'})
vcf.add_info_to_header({'ID': 'INSLEN_ORIG', 'Description': 'Original insertion length', 'Type': 'Integer', 'Number': '.'})
vcf.add_info_to_header({'ID': 'CIPOS95', 'Description': 'Confidence interval (95%) around POS for imprecise variants', 'Type': 'Integer', 'Number': '2'})
vcf.add_info_to_header({'ID': 'CIEND95', 'Description': 'Confidence interval (95%) around END for imprecise variants', 'Type': 'Integer', 'Number': '2'})
vcf.add_info_to_header({'ID': 'SECONDARY', 'Description': 'Secondary breakend in a multi-line variant', 'Type': 'Flag', 'Number': '0'})

max_ins = int(1000)

def convert_variant(v, max_ins):
    set_read_counts(v)
    set_cis_prs(v)
    if v.INFO.get('SVTYPE')=='DEL':
        convert_del(v)
    elif v.INFO.get('SVTYPE')=='DUP':
        convert_dup(v)
    elif v.INFO.get('SVTYPE')=='INV':
        convert_inv(v)
    elif v.INFO.get('SVTYPE')=='INS':
        convert_ins(v, max_ins)
    elif v.INFO.get('SVTYPE')=='BND':
        convert_bnd(v)

def split_ci(ci):
    return[int(ci.split(',')[0]),  int(ci.split(',')[1])]

def uniform_pr(length):
    pr=np.ones(length, dtype='float64')/length
    pr1=','.join( map(str, pr))
    return pr1

def set_read_counts(v):
#    sample=v.sample_list[0]
#    gt=v.genotype(sample)
    pe=0
    sr=0
    if v.format('PR') is not None:# in var.format_dict:
        pe=int(v.format('PR')[0][1])
    if v.format('SR') is not None:# in var.format_dict:
        sr=int(v.format('SR')[0][1])
    v.INFO['PE']=pe
    v.INFO['SR']=sr
    v.INFO['SU']=pe+sr

def set_cis_prs(v):
    imprec=False
    cipos='0,0'
    ciend='0,0'
    prpos=1.0
    prend=1.0
    if v.INFO.get('CIPOS') is not None:# in v.info:
        cipos=v.INFO.get('CIPOS')
        cipos = str(cipos[0]) + ',' + str(cipos[1])
        [start, stop]=split_ci(cipos)
        prpos=uniform_pr(stop-start+1)
        imprec=True
    if v.INFO.get('CIEND') is not None:# in v.info:
        ciend=v.INFO.get('CIEND')
        ciend = str(ciend[0]) + ',' + str(ciend[1])
        [start, stop]=split_ci(ciend)
        prend=uniform_pr(stop-start+1)
        imprec=True
    v.INFO['CIPOS']=cipos
    v.INFO['CIEND']=ciend
    v.INFO['CIPOS95']=cipos
    v.INFO['CIEND95']=ciend
    v.INFO['PRPOS']=prpos
    v.INFO['PREND']=prend
    v.INFO['IMPRECISE'] = imprec

def convert_del(v):
    v.ALT='<DEL>'
    v.INFO['STRANDS']='+-:'+str(v.INFO['SU'])
    v.REF='N'

def convert_dup(v):
    v.ALT='<DUP>'
    v.INFO['STRANDS']='-+:'+str(v.INFO['SU'])
    v.REF='N'

def convert_inv(v):
    v.REF='N'
    strands=''
    if v.INFO.get('INV3') is not None:# in var.info:
        strands='++:'
        v.ALT='N]'+v.CHROM+':'+str(v.INFO['END'])+']'
    else:
        strands='--:'
        v.ALT='['+v.CHROM+':'+str(v.INFO['END'])+'['
    v.INFO['SVTYPE']='BND'
    v.INFO['STRANDS']=strands+str(v.INFO['SU'])    

def convert_ins(v, max_ins):
    v.REF='N'
    v.ALT='<INS>'
    v.INFO['STRANDS']='+.:'+str(v.INFO['SU'])
    orig_len='.'
    new_len=max_ins
    if v.INFO.get('SVLEN') is not None:# in var.info:
        svlen=int(v.INFO.get('SVLEN'))
        orig_len=svlen
        if svlen<max_ins:
            new_len=svlen
    v.INFO['SVLEN']=new_len
    v.INFO['INSLEN_ORIG']=orig_len
        
def convert_bnd(v):
    v.REF='N'
    alt=v.ALT[0]
    ff=alt.find("[")
    newalt=""
    strands=""
    chrom2=""
    breakpoint2=""
    fix_bnd = alt.split(':')
    if '[' in fix_bnd[0]:
        chrom2 = fix_bnd[0].split('[')[1]
    elif ']' in fix_bnd[0]:
        chrom2 = fix_bnd[0].split(']')[1]
    if '[' in fix_bnd[1]:
        breakpoint2 = fix_bnd[1].split('[')[0]
    elif ']' in fix_bnd[1]:
        breakpoint2 = fix_bnd[1].split(']')[0]
    chrom = v.CHROM
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    if chrom2.startswith('chr'):
        chrom2 = chrom2[3:]
    if chrom2 < v.CHROM or (chrom2==v.CHROM and int(breakpoint2)<int(v.POS)):
        v.INFO['SECONDARY'] = True
    if ff==0:
        strands="--:"
        ff1=alt.find("[", 1)
        newalt=alt[0:(ff1+1)]+'N'
    elif ff>0:
        strands="+-:"
        newalt='N'+alt[ff::]
    else:
        ff=alt.find("]")
        if ff==0:
            strands="-+:"
            ff1=alt.find("]", 1)
            newalt=alt[0:(ff1+1)]+'N'
        else:
            strands="++:"
            newalt='N'+alt[ff::]
    v.ALT=newalt
    v.INFO['STRANDS']=strands+str(v.INFO['SU'])

new_vcf = Writer(output_vcf, vcf)
for v in vcf:
    convert_variant(v,max_ins)
    new_vcf.write_record(v)
