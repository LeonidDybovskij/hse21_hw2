import re
from Bio import SeqIO
from collections import defaultdict
from pprint import pprint
from Bio.SeqFeature import SeqFeature, FeatureLocation

mil_1_genome = SeqIO.read("MIL_1.gb", "genbank")

scaffolds = []
for record in SeqIO.parse("scaffolds.fasta", "fasta"):
  scaffolds.append(record)

for i in scaffolds:
  i.annotations['molecule_type'] = 'DNA'
for i in scaffolds:
  i.format('genbank')
SeqIO.write(scaffolds, "GENOME.gbk", "genbank")

with open('gms2.lst', 'r') as gm:
  gm_list = list(gm.readlines())
gm.close()

with open('scaffolds.hits_from_MIL_1.txt', 'r') as mil:
  mil_list = list(mil.readlines())
mil.close()

current_id = 0
mil_dict = defaultdict(str)
for i in mil_list:
  data = i.split('\t')
  if current_id == int(data[0]):
    pass
  else:
    mil_dict[int(data[0])] = data[1][20:]
    current_id = int(data[0])

scaffold_lens = [3533, 7, 3, 1, 1, 1, 3, 1, 1, 5, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 2, 3, 0]
scaffold_lens = scaffold_lens + [1]*10
scaffold_lens = scaffold_lens + [0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0, 0]
scaffold_lens = scaffold_lens + [1]*12
scaffold_lens = scaffold_lens + [0]
scaffold_lens = scaffold_lens + [1]*9

mayfly = 0
sc_dict = {}
for i in gm_list:
  if 'SequenceID' in i:
    if mayfly <=8:
      sc_dict[i[12:28]] = gm_list[(gm_list.index(i) + 1):(gm_list.index(i) + scaffold_lens[mayfly] + 1)]
    else:
      sc_dict[i[12:29]] = gm_list[(gm_list.index(i) + 1):(gm_list.index(i) + scaffold_lens[mayfly] + 1)]
    mayfly += 1

ID = {}
strand = {}
left = {}
right = {}

for key in sc_dict:
  ID[key] = []
  strand[key] = []
  left[key] = []
  right[key] = []
  for i in sc_dict[key]:
    data = i.split()
    ID[key].append(int(data[0]))
    if data[1] == '-':
      strand[key].append(-1)
    else:
      strand[key].append(1)
    if '<' in data[2] or '>' in data[2]:
      left[key].append(int(data[2][1:]))
    else:
      left[key].append(int(data[2]))
    if '<' in data[3] or '>' in data[3]:
      right[key].append(int(data[3][1:]))
    else:
      right[key].append(int(data[3]))

num = 0
for key in sc_dict:
  feature_list = []
  for i in range(len(sc_dict[key])):
    mayfly = SeqFeature(FeatureLocation(left[key][i], right[key][i], strand=strand[key][i]), type="CDS")
    mayfly.qualifiers['locus_tag'] = [str(ID[key][i])]
    mayfly.qualifiers['note'] = mil_dict[ID[key][i]]
    feature_list.append(mayfly)
    # pprint(vars(mayfly))
    # pprint(mayfly.qualifiers)
  for record in SeqIO.parse("proteins.fasta", "fasta"):
    for f in feature_list:
      if record.id == f.qualifiers['locus_tag'][0]:
        f.qualifiers['translation'] = [record.seq]
  for mil_f in mil_1_genome.features:
    if 'protein_id' not in mil_f.qualifiers:
      continue
    for f in feature_list:
      if f.qualifiers['note'] == mil_f.qualifiers['protein_id']:
        f.qualifiers['product'] = mil_f.qualifiers['product']
  scaffolds[num].features = feature_list
SeqIO.write(scaffolds, "GENOME.gbk", "genbank")



'''
scaffold1_cov232 = gm_list[11:3544]
scaffold2_cov196 = gm_list[3548:3555]
scaffold3_cov177 = gm_list[3559:3562]

sc3_f1 = SeqFeature(FeatureLocation(30, 1001, strand=-1), type="CDS")
sc3_f1.qualifiers['locus_tag'] = ['3535']

sc3_f2 = SeqFeature(FeatureLocation(1350, 2192, strand=-1), type="CDS")
sc3_f2.qualifiers['locus_tag'] = ['3536']

sc3_f3 = SeqFeature(FeatureLocation(2304, 2507, strand=1), type="CDS")
sc3_f3.qualifiers['locus_tag'] = ['3537']
pprint(vars(sc3_f1))
record_sc3.features = [sc3_f1, sc3_f2, sc3_f3]
SeqIO.write([record_sc3, record_sc9], "GENOME.gbk", "genbank")

all_feats = [sc3_f1, sc3_f2, sc3_f3]
for record in SeqIO.parse("proteins.fasta", "fasta"):
  for f in all_feats:
    if record.id == f.qualifiers['locus_tag'][0]:
      f.qualifiers['translation'] = [record.seq]

record_sc3.features = [sc3_f1, sc3_f2, sc3_f3]
SeqIO.write([record_sc3, record_sc9], "GENOME.gbk", "genbank")
pprint(vars(sc3_f1))


sc3_f1.qualifiers['note'] = ['CCU72326.1']
sc3_f2.qualifiers['note'] = ['CCU72327.1']
sc3_f3.qualifiers['note'] = ['CCU72328.1']
pprint(sc3_f1.qualifiers)
mil_1_genome = SeqIO.read("/content/T_oleivorans_MIL_1.gbk", "genbank")

for mil_f in mil_1_genome.features:
  if 'protein_id' not in mil_f.qualifiers:
    continue

  for f in all_feats:
    if f.qualifiers['note'][0] == mil_f.qualifiers['protein_id'][0]:
      f.qualifiers['product'] = mil_f.qualifiers['product']

pprint(sc3_f1.qualifiers)
record_sc3.features = [sc3_f1, sc3_f2, sc3_f3]
SeqIO.write([record_sc3, record_sc9], "GENOME.gbk", "genbank")
'''