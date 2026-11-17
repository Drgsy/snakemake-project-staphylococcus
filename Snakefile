# Plik konfiguracyjny - nazwy, foldery etc. 
# W web dev odpowiednik .env, a tutaj envy to odpowiednik package.json w node.js
configfile: 'config.yml'

# Input: genomes_all: 
## Jeden plik to cały genom jednego izolatu bakterii, zapisany jako spakowana FASTA
## Liczba sekwencji w pliku: 1–2, jeśli genom jest kompletny (chromosom i ewentualnie plazmid), albo do kilkuset kontigów, jeśli to draft.
## Skład: 82 genomy S. aureus (różne szczepy) i 18 genomów innych gronkowców, które działają jak kontrola negatywna.

# Input: alleles.fna:
# Warianty sekwencji dla każdego z 7 genów

# Input: bigsdb.tsv:
# Tabela z przypisaniem kombinacji wariantów wszystkich 7 genów do konkretnego szczepu
# Dla każdej kombinacji wariantów wszystkich 7 genów jest przypisany numer szczepu ST

# Wildcards dla nazw plików (w tym przypadku id sekwencji do query i id genomu do subject)
seqids, = glob_wildcards('input/queries/{seqid}.fna')
asmids, = glob_wildcards(config['genomesdir'] + '/{asmid}.fna.gz')

# rule all jest po to by dac znac jakich plików na koncu chcemy, a snakemake sam sobie wyliczy co trzeba zrobić zeby je uzyskać
rule all:
    input:
        'output/merge_mlst/merged_mlst.tsv',
        'output/unique/unique_mlst.fna',
        'output/tree/tree.nwk',
        'output/tree/tree.png'

# Odplanie blasta na kazde querry na kazdy genom (subject) i zapisanie wyników w formacie tsv
# Mamy 100 genomów i 7 sekwencji do wyszukania, więc w sumie 700 (7x100) wyników blasta
rule blastn:
    conda:
        'envs/powb-blast.yml'
    params:
        evalue = 0.01,
        cols   = 'qseqid qlen sseqid length qcovs pident evalue sseq'
    input:
        query = 'input/queries/{seqid}.fna',
        subj  = config['genomesdir'] + '/{asmid}.fna.gz'
    output:
        'output/blastn/blastn_{asmid}_{seqid}.tsv'
    log:
        'log/blastn/blastn_{asmid}_{seqid}.log'
    shell:
        '''gunzip -c {input.subj} 2> {log} |      \
           blastn -evalue         {params.evalue} \
                  -outfmt      "6 {params.cols}"  \
                  -query          {input.query}   \
                  -subject        -               \
                  -out            {output}        \
                   >> {log} 2>&1
        '''

# Ponieważ szukamy w genomie bakteri konkretnych genów spodziewamy się 0-1 odpowiedzi na każde zapytanie.
# Wyjątkiem mogą być przypadki gdzie gen leży na granicy contigów i wtedy może być wykryty w dwóch fragmentach ale następny krok je odrzuci.

# Dla każdego wyniku odpalamy skrypt filter który otrzymaliśmy w zadaniu i filtruje nie interesujące nas wyniki
rule filter:
    conda:
        'envs/powb-python.yml'
    params:
        qcovs = 100,
        cols  = rules.blastn.params.cols
    input:
        rules.blastn.output
    output:
        'output/filter/filtered_{asmid}_{seqid}.tsv'
    log:
        'log/filter/filter_{asmid}_{seqid}.log'
    shell:
        '''scripts/filter.py --eqlen                    \
                             --asmid  {wildcards.asmid} \
                             --qcovs  {params.qcovs}    \
                             --cols  "{params.cols}"    \
                             --input  {input}           \
                             --output {output}          \
                               > {log} 2>&1
        '''

# Łączymy wszystkie przefiltrowane wyniki w jedną tabele tsv dla każdego genomu
# więc otrzymujemy 100 plików wynikowych (dla 100 genomów) każdy z maksymalnie 7 sekwencjami które dopasował blast (głównie 0 lub 7)
rule merge_flt:
    params:
        mask = rules.filter.output[0].replace('{seqid}', '*')
    input:
        expand(rules.filter.output, seqid=seqids, allow_missing=True) # Tutaj tak naprawdę dzieje się grupowanie plików po genomach
    output:
        'output/merge_flt/merged_{asmid}.tsv'
    log:
        'log/merge_flt/merge_{asmid}.log'
    shell:
        # Ten skrypt po prostu łączy pliki w jeden z zachowaniem nagłówka z pierwszego pliku i bez powtarzania nagłówków z kolejnych plików
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail -n +2 ${{file}}
           done >> {output} 2>> {log}
        '''

# Dla każdego zestawu 7 alleli (czyli dla każdego genomu) wyznaczamy typ szczepu na podstawie bazy danych alleli i bazy danych ST
# Skleja sekwencje genów w jeden ciąg
rule mlst:
    conda:
        'envs/powb-python.yml'
    input:
        allelesdb = 'input/alleles.fna',
        mlstdb    = 'input/bigsdb.tsv',
        alleles   = rules.merge_flt.output
    output:
        'output/mlst/mlst_{asmid}.tsv'
    log:
        'log/mlst/mlst_{asmid}.log'
    shell:
        '''scripts/mlst.py --asmid     {wildcards.asmid} \
                           --allelesdb {input.allelesdb} \
                           --mlstdb    {input.mlstdb}    \
                           --alleles   {input.alleles}   \
                           --output    {output}          \
                             > {log} 2>&1
        '''

# Powatrzamy historie z marge_flt tylko tym razem już finalnie łączymy wszystkie wyniki MLST w jeden plik tsv
rule merge_mlst:
    params:
        mask = rules.mlst.output[0].replace('{asmid}', '*')
    input:
        expand(rules.mlst.output, asmid=asmids)
    output:
        'output/merge_mlst/merged_mlst.tsv'
    log:
        'log/merge_mlst/merge_mlst.log'
    shell:
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail -n +2 ${{file}}
           done >> {output} 2>> {log}
        '''

# Uruchamiamy skrypt który z pliku tsv z wynikami MLST:
## Odrzuca genomy z ST = -1 (genomy bez typu szczepu -> niekompletne)
## Usuwamy duplikaty
## Dla każdego szczepu otrzymujemy sekwencje 7 sklejonych genów w jednym pliku fasta
rule unique:
    conda:
        'envs/powb-python.yml'
    input:
        rules.merge_mlst.output
    output:
        'output/unique/unique_mlst.fna'
    log:
        'log/unique/unique.log'
    shell:
        '''scripts/unique.py --input  {input}  \
                             --output {output} \
                               > {log} 2>&1
        '''

# Budowanie drzewa filogenetycznego na podstawie sklejonych sekwencji (odpalanie skryptu)
rule tree:
    conda:
        'envs/powb-python.yml'
    input:
        mlst = rules.merge_mlst.output,
        seqs = rules.unique.output
    output:
        tree = 'output/tree/tree.nwk',
        vis  = 'output/tree/tree.png'
    log:
        'log/tree/tree.log'
    shell:
        '''scripts/tree.py --mlst    {input.mlst}  \
                           --seqs    {input.seqs}  \
                           --outtree {output.tree} \
                           --outvis  {output.vis}  \
                             > {log} 2>&1
        '''

