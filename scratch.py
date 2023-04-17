# found
from Bio import SeqIO
from Bio import Align
from tqdm import tqdm
import pandas as pd
import re



records = SeqIO.parse("./test_files/38872_2#146.ffn", 'fasta')

for record  in records:
    if record.id == 'StNet_00701':
        break

def find_repeats(subject,rep_dictionary):
    spa_type=list()
    for rep in rep_dictionary:
        count = subject.seq.count(rep.seq)
        if(count>0):
            start = 0
            for _ in range(count):
                pos = subject.seq.find(rep.seq,start=start)
                start = pos + len(rep.seq)
                # print(f"{rep.id} found @ {pos}-{start} in {subject.id}")
                spa_type.append({
                    "id":rep.id,
                    "start":pos,
                    "end":start #this could be confusing
                })
    return spa_type

rep_dictionary = list(SeqIO.parse('./database/dictionary.fasta','fasta'))

reps = find_repeats(record,rep_dictionary)

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
hits = []
for rep in tqdm(rep_dictionary):
    alignments = aligner.align("AAAGAAGACAACAAAAACCTGGT", rep.seq)
    for alignment in sorted(alignments):
        # print("Score = %.1f:" % alignment.score)
        hits.append(alignment)

max_score = max([x.score for x in hits])
best = [x for x in hits if x.score==max_score]
len(best)
print(best[0])
type(best[0])
for candidate in best:
    print(candidate)

### build repeat pattern
def region_to_pattern(record,reps,rep_dictionary):
    start = min(reps,key=lambda x:x['start']).get("start")
    end = max(reps,key=lambda x:x['end']).get("end")
    region_with_patterns = str(record.seq[start:end])
    hash_repeats = {}
    for i,rep in enumerate(rep_dictionary):
        hash_repeats[rep.id] = i
    for rep in reps:
        # print(rep.get("id"))
        rep_seq = str(rep_dictionary[hash_repeats[rep.get("id")]].seq)
        region_with_patterns = region_with_patterns.replace(rep_seq,f"-{rep.get('id')}-")
    return region_with_patterns.replace("--","-").strip('-')

pattern = region_to_pattern(record,reps,rep_dictionary)

def suggest_similar_types(pattern,rep_dictionary,database):
    #Check for similar patterns in DB
    # database = pd.read_csv("./database/database.csv",header=0,names=["type","pattern"])
    database = database.dropna()
    #build regex
    regex = re.sub(r"[a-zA-Z]+","",pattern) #erase non numeric
    regex = regex.replace("--","-[^-]*-") #locate missing pattern
    regex = r"" + f"^{regex}$" #build regex
    #extract DNA part
    dna_subject = re.findall(r"-([AGCT]+)-",pattern)
    if len(dna_subject)==0:
        return {}


    df = database[database.pattern.str.match(regex)] #query('pattern.str.contains("r08-r16-r02-r25")')
    #Now get similar types...
    possibles = df.pattern.str.extract(regex)
    hash_repeats = {}
    for i,rep in enumerate(rep_dictionary):
        hash_repeats[rep.id] = i

    #check "distance of subject to target"

    max_score = 0
    hits = {}
    for possible in possibles.iloc[:,0]:
        candidate = rep_dictionary[hash_repeats[f"r{possible}"]]
        alignments = aligner.align(dna_subject[0], candidate.seq)
        for alignament in alignments:
            # print(f"New best: {alignament.score}")
            if alignament.score >= max_score:
                max_score = alignament.score
                hits[f"{possible}"] = alignament.score
    suggestions = {}
    for k in hits:
        suggest = f"08-16-02-25-{k}-25"
        suggestions[database.query("pattern==@suggest")]=suggest
    
    return suggestions