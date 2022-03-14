import vcf2excel
import train_model

if __name__ == '__main__':

    ## 循环运行文件
    # for i in range(1,21):
    #     uni_dir='D:/2_Professional Courses/132_基因项目/GWAS数据/新增_20基因数据/'
    #     file_address=uni_dir+str(i)+'.vcf'
    #     print("目前正在运行文件：  "+ str(i)+'.vcf')
    #     train_model.train_model(file_address)

    ## 单独运行文件
    user_file=input('请输入用户文件路径： \n')
    train_model.train_model(user_file)


(perl -alne '{print if /^rs/}' dm_23andme_v3_110219.txt  |cut -f 1 >23andme.rsID.listcat ../variation/autochr.highQuali.dbsnp.vcf  23andme.rsID.list |perl -alne '{if($F[2]=~/^rs/){if(/1\/1/){$gt=$F[4].$F[4]}else{$gt=$F[3].$F[4]};$h{$F[2]}="$F[0]\t$F[1]\t$gt" }  print "$_\t$h{$_}" if /^rs/}' >my_23andme.1.txtzcat ~/annotation/variation/human/dbSNP/All_20160601.vcf.gz |perl -alne 'BEGIN{ open FH,"my_23andme.1.txt";while(<FH>){chomp;@F=split;if(/^rs/){ $pos{$.}=$_;if($F[3]){$h{$F[0]}=$_}else{$tmp{$F[0]}=1}  }} }{if(exists $tmp{$F[2]}){ $tmp{$F[2]}="$F[0]\t$F[1]\t$F[2]$F[2]"  }}END{foreach(sort{$a<=>$b} keys %pos){ if(exists $h{$pos{$_}} ){$value=$h{$pos{$_}}}else{$value=$tmp{$pos{$_}} } ;print "$pos{$_}\t$value" }}'







#model.SC_model('D:/2_Professional Courses/132_基因项目/GWAS数据/新增_20基因数据/1.vcf')


