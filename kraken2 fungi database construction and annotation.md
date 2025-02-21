# 利用kraken2构建真菌数据库，真菌参考基因组来自于NCBI RefSeq(<https://www.ncbi.nlm.nih.gov/refseq>), FungiDB(<http://fungidb.org>), Ensemble(<http://fungi.ensembl.org>)和最新发表的cultivated gut fungi(GCFs) (<https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833221>)，利用kraken2构建自定义数据库的细节可参考<https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown>

## 1.下载fungi参考基因组

### 1.1 NCBI RefSeq fungi:  <https://www.ncbi.nlm.nih.gov/refseq>

在构建真菌数据库时可直接下载kraken2已经构建好的数据库中的fungi部分，如果要最新版本的fungi数据，可参考这个代码进行下载

```bash
cd /fungi/fungi_ncbi

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt -O assembly_summary_fungi.txt

wc -l assembly_summary_fungi.txt #634,带两行表头和 行无效链接，实际 

awk -F"\t" '{print $20}' assembly_summary_fungi.txt > fungi_ftplinks.txt
#查看有没有无效的下载链接
grep "na_genomic.fna.gz" fungi_ftplinks.txt  | wc -l #0

# Append fna suffix to get full download link，为ftp路径添加后缀
awk 'BEGIN{FS=OFS="/"; filesuffix="genomic.fna.gz"} NR==1 {print} NR>1 {ftpdir=$0; asm=$NF; file=asm"_"filesuffix; print ftpdir"/"file}' fungi_ftplinks.txt > tmpfile.txt && mv tmpfile.txt fungi_ftplinks.txt

######## 并行下载真菌基因组文件
cd ~/data5/baidefeng/fungi/fungi_ncbi
#parallel -j 50 wget < ../bacterial_pre2020_ftplinks.txt
#确保文件夹内文件数为187974，如果不是，请重复运行此命令
cat ~/fungi/fungi_ncbi/fungi_ftplinks.txt | xargs -n 1 -P 50 wget -c
##查看文件数
ls -1 | wc -l  #列出当前文件夹内的所有文件和子文件夹数
find . -type f | wc -l #只看文件数187974

# 解压缩fna.gz文件
cd ~/fungi/fungi_ncbi
gunzip *.fna.gz
```

### 1.2 FungiDB(<http://fungidb.org>)

FungiDB数据库中的真菌基因组需要到<http://fungidb.org>进行下载

### 1.3 Ensemble(<http://fungi.ensembl.org>)

Ensemble数据库中的真菌基因组需要到<http://fungi.ensembl.org>进行下载

### 1.4 Cultivated gut fungi (GCFs) (<https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833221>)

```bash
# 支持断点续传利用项目号从NCBI上下载参考基因组
# 获取项目的数据集 
datasets download genome accession PRJNA833221 --dehydrated --include genome 

# 解压下载的文件 
unzip ncbi_dataset.zip -d ncbi_dataset cd ncbi_dataset 

# 下载实际数据文件：使用 datasets rehydrate 命令下载所有相关的 FASTA 文件： 
# 支持断点续传，下载中断了请重复运行 
datasets rehydrate --directory .
```

## 2. 参考基因组格式标准化

### 2.1 根据物种名在NCBI上批量提取TaxID和Taxonomy信息

```python
from Bio import Entrez
import csv
import pandas as pd

# 设置邮箱（NCBI要求提供）
Entrez.email = "baidefeng1234win@163.com"

# 从NCBI检索物种的分类信息
def get_taxonomy_by_name(species_name):
    # 使用esearch查找物种的TaxID
    search_handle = Entrez.esearch(db="taxonomy", term=species_name)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    # 如果找不到匹配的物种，返回None
    if not search_results["IdList"]:
        return None
    
    taxid = search_results["IdList"][0]
    
    # 使用efetch来获取分类信息
    fetch_handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(fetch_handle)
    fetch_handle.close()
    
    # 提取分类信息和物种名称
    record = records[0]
    lineage = record["LineageEx"]
    
    # 构造分类信息字符串
    classification = []
    for taxon in lineage:
        rank = taxon["Rank"]
        name = taxon["ScientificName"]
        if rank == "superkingdom":
            classification.append(f"d__{name}")
        elif rank == "kingdom":
            classification.append(f"k__{name}")
        elif rank == "phylum":
            classification.append(f"p__{name}")
        elif rank == "class":
            classification.append(f"c__{name}")
        elif rank == "order":
            classification.append(f"o__{name}")
        elif rank == "family":
            classification.append(f"f__{name}")
        elif rank == "genus":
            classification.append(f"g__{name}")
        elif rank == "species":
            classification.append(f"s__{name}")
    
    # 加入最后一级的物种具体名字
    species_specific = f"s1__{record['ScientificName'].replace(' ', '_')}"
    classification.append(species_specific)
    
    return taxid, ";".join(classification)

# 读取输入CSV文件，提取物种名称列
input_csv = "taxonomy_info.csv"  # 输入的CSV文件
output_csv = "taxonomy_info_out4.csv"  # 输出的CSV文件
species_column = "Species Name"  # 物种名称列

# 使用pandas读取输入的CSV文件
df = pd.read_csv(input_csv)

# 打开输出文件并写入
with open(output_csv, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Species Name", "TaxID", "Classification"])  # 写入表头
    
    # 遍历物种名称列表
    for species in df[species_column]:
        result = get_taxonomy_by_name(species)
        if result:
            taxid, classification = result
            writer.writerow([species, taxid, classification])
        else:
            writer.writerow([species, "Not Found", ""])

print(f"分类信息已保存到 {output_csv}")
```

### 2.2 修改下载得到的fasta/fa/fna文件的head，使其符合kraken2构建数据库的要求

```bash
# 单个修改fna文件的head，如果是fa或者fasta文件，后缀名修改成对应的格式
awk '/^>/ { sub(">", "", $1); $0 = ">" $1 "|kraken:taxid|5065 Aspergillus sp. S4402 ctg_1" } 1' GCA_023619285.1_ASM2361928v1_genomic.fna>GCA_023619285.1_ASM2361928v1_genomic_modified.fna

# 批量修改fna文件的head
# 准备一个名称为fna_info.tsv文件的文件,包含三列信息：filename, taxid, organism_name, filename列是文件名称，taxid是NCBI上对应物种的NCBI Taxonomy ID, organism_name是物种名称
# 跳过CSV文件的标题行，遍历数据表中的每一行
# fna文件，如果是fa或者fasta文件，后缀名修改成对应的格式
tail -n +2 Cell_fna_info.csv | while IFS=',' read -r filename taxid organism_name; do
    # 如果文件存在
    if [[ -f "$filename" ]]; then
        # 使用awk修改文件头信息
        awk -v taxid="$taxid" -v organism="$organism_name" '/^>/ {
            sub(">", "", $1);
            $0 = ">" $1 "|kraken:taxid|" taxid " " organism
        } 1' "$filename" > "${filename%.fna}_modified.fna"
        
        echo "Modified $filename and saved as ${filename%.fna}_modified.fna"
    else
        echo "File $filename not found!"
    fi
done

# 整合数据
# 将选择的fna文件合并
# Cultivated gut fungi (GCFs)
cat *modified.fna > Cell_combined.fna

# FungiDB数据库
cat *modified.fasta > FungiDB_combined.fna

# Ensemble数据库
cat *modified.fa > Ensemble_combined.fna

# Cultivated gut fungi (GCFs); FungiDB; Ensemble三者合并
cat /data5/baidefeng/fungi/FungiDB2/FungiDB_combined.fna /data5/baidefeng/fungi/Ensemble/Ensemble_combined.fna /data5/baidefeng/fungi/ncbi_dataset/Cell_combined.fna > /data5/baidefeng/fungi/FungiDB_ensemble_cell_merged01.fna
```

## 3.利用kraken2构建fungi数据库

```bash
# 利用conda安装kraken2和bracken
conda create -n kraken2-2.1.3
conda activate kraken2-2.1.3
conda install -c bioconda kraken2=2.1.3
conda install -c bioconda bracken

# 进入kraken2环境
conda activate kraken2-2.1.3

# 下载kraken2数据库中的fungi部分
kraken2-build --download-library fungi --db kraken2db2 --threads 12
kraken2-build --download-taxonomy --db ./kraken2db2 --threads 24 

# 建库测试
kraken2-build --build --db ./kraken2db --threads 36

# 将fna文件添加到数据库中
conda activate kraken2-2.1.3
kraken2-build --add-to-library /data5/baidefeng/fungi/FungiDB_ensemble_cell_merged01.fna --db ./kraken2db --threads 36

# 建库
kraken2-build --build --db ./kraken2db --threads 36
```

## 4.真菌物种注释

```bash
# 单个注释, 输入的文件是质控和去宿主后的fastq文件
# 这个也运行正常，奇怪
kraken2 \
  --db /data5/baidefeng/kraken2db2 \
  --threads 12 --use-names --report-zero-counts \
  --paired /data5/baidefeng/npc/temp/hr/R056_paired_1.fastq /data5/baidefeng/npc/temp/hr/R056_paired_2.fastq \
  --report /data5/baidefeng/npc/temp/kraken2_fungi3/R056_3.report \
  --output /data5/baidefeng/npc/temp/kraken2_fungi3/R056_3.output

python3 kreport2mpa.py -r /data5/baidefeng/npc/temp/kraken2_fungi3/R056_3.report \
    --display-header -o /data5/baidefeng/npc/temp/kraken2_fungi3/R056_3.mpa

mkdir -p age
cd age
mkdir -p seq temp result

cp /data1/liuyongxin/age/meta/result/metadata.txt /data5/baidefeng/age/result/

# 批量注释, 输入的文件是质控和去宿主后的fastq文件
# 元数据细节优化
# 转换Windows回车为Linux换行
sed -i 's/\r//' result/metadata01.txt

# 去除数据中的一个多余空格
sed -i 's/Male  /Male/' result/metadata01.txt

# 百岁老人
#  多样本并行生成report，1样本8线程逐个运行，内存大但速度快，不建议用多任务并行
for i in `tail -n+2 /data5/baidefeng/age/result/metadata_sardinian.txt | cut -f1`;do
  kraken2 --db /data5/baidefeng/kraken2db2 \
  --paired /data3/liyahui/age/sardinian/temp/hr/${i}_1.fastq /data3/liyahui/age/sardinian/temp/hr/${i}_2.fastq \
  --threads 8 --use-names --report-zero-counts \
  --report /data5/baidefeng/age/temp/kraken2_fungi_sardinian/${i}.report \
  --output /data5/baidefeng/age/temp/kraken2_fungi_sardinian/${i}.output; done

# 使用krakentools转换report为mpa格式
for i in `tail -n+2 /data5/baidefeng/age/result/metadata_sardinian.txt | cut -f1`;do
  python3 kreport2mpa.py -r /data5/baidefeng/age/temp/kraken2_fungi_sardinian/${i}.report \
    --display-header -o /data5/baidefeng/age/temp/kraken2_fungi_sardinian/${i}.mpa; done

# 合并样本为表格
mkdir -p result/kraken2_fungi_sardinian_test1

# 输出结果行数相同，但不一定顺序一致，要重新排序
tail -n+2 /data5/baidefeng/age/result/metadata_sardinian.txt | cut -f1 | rush -j 1 \
  'tail -n+2 /data5/baidefeng/age/temp/kraken2_fungi_sardinian/{1}.mpa | LC_ALL=C sort | cut -f 2 | sed "1 s/^/{1}\n/" > /data5/baidefeng/age/temp/kraken2_fungi_sardinian/{1}_count '

# 提取第一样本品行名为表行名
header=`tail -n 1 /data5/baidefeng/age/result/metadata_sardinian.txt | cut -f 1`
echo $header
tail -n+2 /data5/baidefeng/age/temp/kraken2_fungi_sardinian/${header}.mpa | LC_ALL=C sort | cut -f 1 | \
  sed "1 s/^/Taxonomy\n/" > /data5/baidefeng/age/temp/kraken2_fungi_sardinian/0header_count
head -n3 /data5/baidefeng/age/temp/kraken2_fungi_sardinian/0header_count

# paste合并样本为表格
ls /data5/baidefeng/age/temp/kraken2_fungi_sardinian/*count
paste /data5/baidefeng/age/temp/kraken2_fungi_sardinian/*count > /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/tax_count.mpa

# 检查表格及统计
csvtk -t stat /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/tax_count.mpa
head -n 5 /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/tax_count.mpa


# bracken安装
conda install -y bracken
bracken -h

# Generate the Bracken database file (databaseXmers.kmer_distrib)
# 遇到bug:/data5/baidefeng/miniconda3/bin/bracken-build: line 231: syntax error: unexpected end of file
# 解决方法：在https://github.com/jenniferlu717/Bracken下载bracken-build对/data5/baidefeng/miniconda3/bin/bracken-build进行替换
bracken-build -d /data5/baidefeng/kraken2db2 -t 16 -k 35 -l 150

# Bracken丰度估计
mkdir -p /data5/baidefeng/age/temp/bracken_fungi_sardinian

conda activate kraken2.1.3

# 测序数据长度，通常为150，早期有100/75/50/25
readLen=150

# 20%样本中存在才保留
prop=0.05

# 设置分类级D,P,C,O,F,G,S，常用界D门P和属G种S
for tax in D P G S;do
tax=S
for i in `tail -n+2 /data5/baidefeng/age/result/metadata.txt | cut -f1`;do
    # i=C1
    bracken -d /data5/baidefeng/kraken2db2 \
      -i /data5/baidefeng/age/temp/kraken2_fungi1/${i}.report \
      -r ${readLen} -l ${tax} -t 12 \
      -o /data5/baidefeng/age/temp/bracken_fungi1/${i}.brk \
      -w /data5/baidefeng/age/temp/bracken_fungi1/${i}.report; done

tax=S
for i in `tail -n+2 /data5/baidefeng/age/result/metadata_sardinian.txt | cut -f1`;do
    # i=C1
    bracken -d /data5/baidefeng/kraken2db2 \
      -i /data5/baidefeng/age/temp/kraken2_fungi_sardinian/${i}.report \
      -r ${readLen} -l ${tax} -t 12 \
      -o /data5/baidefeng/age/temp/bracken_fungi_sardinian/${i}.brk \
      -w /data5/baidefeng/age/temp/bracken_fungi_sardinian/${i}.report; done


# 需要确认行数一致才能按以下方法合并      
#wc -l /data5/baidefeng/age/temp/bracken_fungi1/*.report
wc -l /data5/baidefeng/age/temp/bracken_fungi_sardinian/*.report

# 行名不一致
# 提取所有report文件中独特的物种名
# 假设所有的 .report 文件都在当前目录下
for file in /data5/baidefeng/age/temp/bracken_fungi_sardinian/*.brk; do cut -f1 $file; done | sort | uniq > /data5/baidefeng/age/result/species_list_sardinian.txt

# 统一每个brk文件的行数
for file in /data5/baidefeng/age/temp/bracken_fungi_sardinian/*.brk; do
    awk 'FNR==NR {species[$1]; next} {if($1 in species) print $0; else print $1, "0", "0", "0", "0", "0"}' /data5/baidefeng/age/result/species_list_sardinian.txt $file > temp_file && mv temp_file $file
done


# 提取所有 .brk 文件中的唯一物种名，生成统一的 species_list.txt
for file in /data5/baidefeng/age/temp/bracken_fungi_sardinian2/*.brk; do
    cut -f1 "$file"
done | sort | uniq > /data5/baidefeng/age/result/species_list_sardinian.txt

# 确保目标目录存在
mkdir -p temp/bracken_fungi_sardinian2/

# 批量修改
# 定义文件路径
SPECIES_FILE="/data5/baidefeng/age/result/species_list_sardinian.txt"
OUTPUT_DIR="/data5/baidefeng/age/temp/bracken_fungi_sardinian2"

# 遍历所有的 .brk 文件
for EB55_FILE in temp/bracken_fungi_sardinian2/*.brk; do
    # 获取文件名，不包括路径
    FILENAME=$(basename "$EB55_FILE")
    OUTPUT_FILE="$OUTPUT_DIR/${FILENAME%.brk}_updated.brk"

    # 提取现有物种名
    cut -f1 "$EB55_FILE" > "$OUTPUT_DIR/existing_species.txt"

    # 查找缺失的物种名
    grep -Fxv -f "$OUTPUT_DIR/existing_species.txt" "$SPECIES_FILE" > "$OUTPUT_DIR/missing_species.txt"

    # 创建填充0的行，格式为：物种名后跟五个0
    awk -F'\t' '{print $1 "\t0\t0\t0\t0\t0"}' "$OUTPUT_DIR/missing_species.txt" > "$OUTPUT_DIR/new_lines.txt"

    # 将原始的 .brk 文件和新增的行合并，并生成最终文件
    cat "$EB55_FILE" "$OUTPUT_DIR/new_lines.txt" > "$OUTPUT_FILE"

    # 清理中间文件
    rm "$OUTPUT_DIR/existing_species.txt" "$OUTPUT_DIR/missing_species.txt" "$OUTPUT_DIR/new_lines.txt"
done

# bracken结果合并成表: 需按表头排序，提取第6列reads count，并添加样本名
tail -n+2 /data5/baidefeng/age/result/metadata_sardinian.txt | cut -f1 | rush -j 1 \
  'tail -n+2 /data5/baidefeng/age/temp/bracken_fungi_sardinian2/{1}_updated.brk | LC_ALL=C sort | cut -f6 | sed "1 s/^/{1}\n/" \
  > /data5/baidefeng/age/temp/bracken_fungi_sardinian2/{1}.count'

# 提取第一样本品行名为表行名
h=`tail -n1 /data5/baidefeng/age/result/metadata_sardinian.txt|cut -f1`
tail -n+2 /data5/baidefeng/age/temp/bracken_fungi_sardinian2/${h}_updated.brk | LC_ALL=C sort | cut -f1 | \
  sed "1 s/^/Taxonomy\n/" > /data5/baidefeng/age/temp/bracken_fungi_sardinian2/0header.count

# 检查文件数，为n+1
ls /data5/baidefeng/age/temp/bracken_fungi_sardinian2/*count | wc

# paste合并样本为表格，并删除非零行
paste /data5/baidefeng/age/temp/bracken_fungi_sardinian2/*count > /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.${tax}_new01.txt

# 统计行列，默认去除表头
csvtk -t stat /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.${tax}_new01.txt

db=/data/meta/db

# 按频率过滤，-r可标准化，-e过滤(microbiome_helper)
Rscript ${db}/EasyMicrobiome/script/filter_feature_table.R \
  -i result/kraken2_fungi_all_zhang2022/bracken.${tax}_new01.txt \
  -p ${prop} \
  -o result/kraken2_fungi_all_zhang2022/bracken.${tax}.${prop}
# head result/kraken2/bracken.${tax}.${prop}
# done

Rscript ${db}/EasyMicrobiome/script/filter_feature_table.R \
  -i /data5/baidefeng/age/result/kraken2_origin_self2/bracken.${tax}_new01.txt \
  -p ${prop} \
  -o /data5/baidefeng/age/result/kraken2_origin_self2/bracken.${tax}.${prop}
# head result/kraken2/bracken.${tax}.${prop}
# done

csvtk -t stat /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.?_new01.txt
csvtk -t stat /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.?.$prop

# 个性化结果筛选
# 门水平去除脊索动物(人)
grep 'Chordata' /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.S.${prop}
grep -v 'Chordata' /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.S.${prop} > /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.S.${prop}-H

# 按物种名手动去除宿主污染，以人为例(需按种水平计算相关结果)
# 种水平去除人类P:Chordata,S:Homo sapiens
grep 'Homo sapiens' /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.S.${prop}
grep -v 'Homo sapiens' /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.S.${prop} \
  > /data5/baidefeng/age/result/kraken2_fungi_sardinian_test1/bracken.S.${prop}-H

# 分析后清理每条序列的注释大文件
/bin/rm -rf /data5/baidefeng/age/temp/kraken2_fungi_sardinian/*.output
```

