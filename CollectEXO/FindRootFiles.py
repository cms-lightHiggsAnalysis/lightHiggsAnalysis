import subprocess
import sys
import fileinput
import os
import re
def getFileList(inDir, inputName):
    tempfile = open('tmpFile.txt','w')
    cmd_all='find ' + inDir + ' -maxdepth 2 -mindepth 2 -name *Version1' 
    print cmd_all
    process_all=subprocess.Popen(cmd_all.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    tempfile.write(process_all.communicate()[0])
    tempfile.close()
    outDirList=[]
    with open('tmpFile.txt','r') as readfile:
        outDirList=[line.rstrip() for line in readfile]
    readfile.close()
    os.remove('tmpFile.txt')
    JobList=[]
    for Dir in outDirList:
        tempfile_Dir=open('tmpDirFile.txt','w')
        cmd_last='find ' + Dir +' -mindepth 2 -maxdepth 2'
        process_last=subprocess.Popen(cmd_last.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        tempfile_Dir.write(process_last.communicate()[0])
        tempfile_Dir.close()
        last_line=subprocess.check_output(['tail','-1', 'tmpDirFile.txt']).rstrip()
        JobList.append(last_line)
        os.remove('tmpDirFile.txt')
    for Job in JobList:
        print Job
        result=Job.split('/')[8]#6 when inDir_ is /eos/cms/store/user/mshi
        fileName='./AllRootFiles/'+result+'.txt'
        getFileInOneFolderList(Job, fileName)
        
def getFileInOneFolderList(inDir, inputName):
    tempfile = open(inputName,'w')
    cmd = 'find ' + inDir + ' -name *TFile*.root'
    print cmd
    process = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    tempfile.write(process.communicate()[0])
    tempfile.close()
    for line in fileinput.input([inputName], inplace=True):
        sys.stdout.write('root://eoscms/{l}'.format(l=line))


inputName_='testscript.txt'
inDir_='/eos/cms/store/group/phys_higgs/HiggsExo/mshi'
#inDir_= '/eos/cms/store/user/mshi'
getFileList(inDir_,inputName_)

