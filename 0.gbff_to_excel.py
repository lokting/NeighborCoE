#提取gbff文件信息为excel文件

# 导入 SeqIO 模块和 FeatureIO 模块
import os
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import openpyxl as xl

# 读取 GBFF 文件
# 输入路径
gbff_dir = r"input/genome_gbff/"
gbff_files = [f for f in os.listdir(gbff_dir) if f.endswith(".gbff")]
gbff_file = os.path.join(gbff_dir, gbff_files[0])
# 输出路径
excel_path = r"tmp_output/gbff_function.xlsx"

records = SeqIO.parse(gbff_file, "genbank")

# 创建 Excel 表格
workbook = xl.Workbook()
sheet = workbook.active
sheet.title = "Sheet1"
sheet.append(["locus_tag", "product","translation",  "protein_id", "gene", "location_start","location_end","Strand","pseudo"])

# 遍历 SeqRecord 对象
for record in records:
    # 遍历 features 属性
    for feature in record.features:
        # 判断 feature 是否为 CDS 类型
        if feature.type == "CDS":

            information_list = []

            # 获取 locus_tag、translation 等信息
            gene = feature.qualifiers.get("gene", [""])[0]
            product = feature.qualifiers.get("product", [""])[0]
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            protein_id = feature.qualifiers.get("protein_id", [""])[0]
            translation = feature.qualifiers.get("translation", [""])[0]
            pseudo = feature.qualifiers.get("pseudo", "0")[0]

            # 获取 CDS 的位置信息
            location = feature.location
            # start = location.start.position + 1
            start = location.start + 1
            # end = location.end.position
            end = location.end
            
            if location.strand == 1 :
                strand = "+"
            elif location.strand == -1:
                strand = "-"

            if pseudo == "":
                pseudo_output = True
            elif pseudo == "0":
                pseudo_output = False
            # 输出 CDS 的相关信息
            # print("location:", location)
            # print("gene:", gene)
            # print("locus_tag:", locus_tag)
            # print("product:", product)
            # print("protein_id:", protein_id)
            # print("translation:", translation)

            #添加信息
            information_list.append(locus_tag)
            information_list.append(product)
            information_list.append(translation)
            information_list.append(protein_id)
            information_list.append(gene)
            information_list.append(str(start))
            information_list.append(str(end))
            information_list.append(str(strand))
            information_list.append(str(pseudo_output))
            sheet.append(information_list)

workbook.save(excel_path)

