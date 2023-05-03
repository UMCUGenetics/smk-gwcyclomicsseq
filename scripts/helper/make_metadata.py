# This script convert csv copied from google drive 20230321-GWC-metadata-Dami_v1.0 "RunName" and "SampleName" columns and convert them to desired samplesheet.csv format.

# Select a genome reference

# genome_reference = "GCh38.p14"
genome_reference = "hs37d5"

# Filename exported from google sheet
file_in = "/Users/liting/Desktop/hg37.csv"

# outfile
file_out = "/Users/liting/01_data/20230407/samplesheet_hg37.csv"

# Parsing
with open(file_in) as f:
    b = f.readlines()
    
b = [line.split(",") for line in b]
filenames = []
run_names = []
sample_names = []
for x in range(len(b)):
    run_name = b[x][0].strip("\"") 
    sample_name = b[x][1].strip()
    if run_name == 'DER3846':
        continue
    print(run_name)
    filenames.append( sample_name +  "_"+ run_name + ".YM_gt_3.bam")
    sample_names.append(sample_name)
    run_names.append(run_name)

df = pd.DataFrame([filenames, sample_names, run_names]).T
df.columns = ['bam_name', 'sample_name', 'run_name' ]
df['genome_reference'] = genome_reference
df = df.sort_values('bam_name')
# Check if data is correct

# Write CSV
df.to_csv(file_out, index=False)
