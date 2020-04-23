import sys, os 
import chimera
from chimera import openModels
from chimera import runCommand

prots = ['CSN1', 'CSN2', 'CSN3', 'CSN4',
         'CSN5', 'CSN6', 'CSN7', 'CSN8']

chains = {'CSN1':"A",\
          'CSN2':"B",\
          'CSN3':"C",\
          'CSN4':"D",\
          'CSN5':"E",\
          'CSN6':"F",\
          'CSN7':"G",\
          'CSN8':"H"
}

vol = {'CSN1':67195,\
       'CSN2':62428,\
       'CSN3':25520,\
       'CSN4':55981,\
       'CSN5':45468,\
       'CSN6':50423,\
       'CSN7':35839,\
       'CSN8':28101
       }

col =  {'CSN1':"#800080",\
        'CSN2':"#0000ff",\
        'CSN3':"#ff3333",\
        'CSN4':"#6666ff",\
        'CSN5':"#474747",\
        'CSN6':"#008c9c",\
        'CSN7':"#ffd900",\
        'CSN8':"#228c22"
        }

totvol = 0.0
for p in prots:
    totvol += vol[p]

suffix='DSSO_DHSO_BMSO'
os.chdir("./Structure_%s" % suffix)
print(suffix)

#from chimera.tkgui import saveReplyLog
#saveReplyLog("./reply-log.txt")

runCommand('set bgcolor white')

i=0
#runCommand('open cluster_center_model.rmf3')
#i += 1
# Whole LPD
runCommand('open LPD_Whole.mrc')
runCommand('color #000000 #' +str(i))
runCommand('volume #'+str(i)+' step 1 transparency 0.70 encloseVolume ' + str(1.0*totvol))
i+=1    

# Read components
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' encloseVolume '+str(1.5*vol[p])+' step 1 transparency 0.70') 
    i += 1
    
runCommand('open ./%s.pdb' % suffix)
for c in prots:
    runCommand('color '+ col[c] + ' #'+str(i)+':.'+chains[c])

runCommand('molmap #'+str(i)+' 30')


#runCommand('open CSN.pdb')
#runCommand('del #'+str(i)+':.I-Z')

runCommand('vop add #1-8')
runCommand('fitmap #10 #9.1 search 20 listFits false')

runCommand('matrixcopy #10 #9.1 moving #0-8')

runCommand('~modeldisplay #0 #9.1 #10')
runCommand('modeldisplay #1-8')

runCommand('center')

runCommand('vop resample #0 onGrid #9.1')
runCommand('volume #11 save '+suffix+'.Whole.LPD.mrc')
#runCommand('matrixcopy 9 10')
for i in range(1,9):
  runCommand('vop resample #'+str(i)+' onGrid #9.1')
  runCommand('volume #'+str(11+i)+' save '+ suffix+'.'+prots[i-1]+'.LPD.mrc')

#runCommand('vop resample #10 onGrid #9.1')
#runCommand('volume #10 save '+ suffix+'.Added.LPD.mrc')

exit()
