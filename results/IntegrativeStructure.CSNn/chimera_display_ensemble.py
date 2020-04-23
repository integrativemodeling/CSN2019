import chimera
from chimera import openModels
from chimera import runCommand
import sys, os 
prots = ['CSN1', 'CSN2', 'CSN3', 'CSN4',
         'CSN5', 'CSN6', 'CSN7', 'CSN8', 'CSN9',
         'CSN9120', 'CSN92040', 'CSN94057']

prots0 = ['CSN1', 'CSN2', 'CSN3', 'CSN4',
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
           'CSN8':28101,\
           'CSN9':7515,\
           'CSN9120':2505,\
           'CSN92040':2505,\
           'CSN94057':2505
       }

col =  {'CSN1':"#800080",\
            'CSN2':"#0000ff",\
            'CSN3':"#ff3333",\
            'CSN4':"#6666ff",\
            'CSN5':"#474747",\
            'CSN6':"#008c9c",\
            'CSN7':"#ffd900",\
            'CSN8':"#228c22", \
            'CSN9':"#FF8C00", \
            'CSN9120':"#FF7B00",\
            'CSN92040':"#FF9000",\
            'CSN94057':"#FFA500"
        }

precision = 2.2

totvol = 0.0
for p in prots:
    totvol += vol[p]

suffix='DSSO_DHSO_BMSO'
os.chdir("./Structure_%s" % suffix)

print(suffix)
runCommand('set bgcolor white')
i=0
#runCommand('open cluster_center_model.rmf3')
#i += 1
# Whole LPD
runCommand('open '+suffix+'.Whole.LPD.mrc')
runCommand('color #000000 #'+str(i))
runCommand('volume #'+str(i)+' step 1 transparency 0.70 encloseVolume ' + str(precision * totvol))
i+=1    


# Read components
for p in prots:
    runCommand('open '+suffix+'.'+p+'.LPD.mrc')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' encloseVolume '+str(precision*vol[p])+' step 1 transparency 0.70')
    runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' style bold size 52 xpos .1 ypos ' + str( 0.7 - i / 20.0)) 
    i += 1
    
runCommand('open ./%s.pdb' % suffix)
print(i)
#runCommand('open CSN.pdb')
#runCommand('del #'+str(i)+':.I-Z')

for c in prots0:
    runCommand('color '+col[c]+' #'+str(i)+':.'+chains[c])


runCommand('~modeldisplay #0')

