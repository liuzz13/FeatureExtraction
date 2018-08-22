#输入 待转换textgrid的wavlist，在textgrid所在文件夹下生成提取后的文字 .txt

import sys
import re

#f1 = open(sys.argv[1],"r")
f1 = open("E:\\培训\\语言模型\\wavlist.txt","r")
filelist = f1.readlines()

m = "哎不客气，那稍后麻烦做评价，感谢来电，女士再见\n"
pattern = re.compile("\" [\u4e00-\u9fa5]",re.U)
n = re.match(pattern,m) != None
print (n)


for file in filelist:
	f = open(file[:-1],"r",encoding='utf-16-be')
	f2 = open(file[:-8] + ".txt","w")
	text = f.readlines()
	for line in text:
		#print (line)
		if (re.match('\"[\u4e00-\u9fa5]',line) != None):	#利用正则表达式匹配
			f2.write(line)

