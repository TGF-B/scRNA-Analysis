#ALWAYS remember harware disk(/media/...) is for data storage, 
#while local disk(/home/...) is for software installation,
#so you can avoid many unnecessary link parsing problems

#1:install the linux(ubuntu 20.04 ) with rufus, download rufus and ubuntu.iso locally, 

#2:write the ubuntu.iso into usb disk in FAT32 format(if not, re-format it) by rufus.

#3:plug in the usb, restart the computer, press F7(might vary) to enter the boot menu, select your usb disk to boot.

#4:waiting for ubuntu.iso to  be installed.(1 hour at least)

#5:Once Ubuntu desktop lauched, intall softwares whichever you want.

#6:install miniconda(suggest) locally, NOT in hard disk!

#7:install cellranger locally, NOT in hard disk!

#8:set up the PATH environment for cellranger, so you can launch cellranger by simply cellranger command.

#9:run the cellranger.sh script to process the fastq files.
#(pay attention to the file paths and functions you need, which you could refer to 10X official instructions)

#other tips:
#1:Cellranger 9.0.0 support automatic cell annotation which by Cloud(CLI), but calling for 10X Chromium token,
#somtimes could encounter network problems( I failed every time even with), so you better use the previous version of cellranger(eg. 8.0.0/8.0.1)

#2:Download specific version of cellranger by official command.

#3:Usually we put the cellranger count command in a shell script, using chmod +x xxx.sh and ./xxx.sh,
#sometimes it would return permission denied errors, under such circumstances, you can use sudo ./xxx.sh

#cellranger count
chmod +x /media/donaldtangai4s/Yu_Omics2/fastq/cellranger.sh
#if you do the second run 
sudo rm /home/donaldtangai4s/new/Ctrl-3/_lock 

sudo /media/donaldtangai4s/Yu_Omics2/fastq/cellranger.sh


#cellranger annotate

chmod +x /media/donaldtangai4s/Yu_Omics2/fastq/annotate.sh

sudo /media/donaldtangai4s/Yu_Omics2/fastq/annotate.sh

#scVelo
chmod +x /media/donaldtangai4s/Yu_Omics2/fastq/loom.sh

sudo /media/donaldtangai4s/Yu_Omics2/fastq/loom.sh