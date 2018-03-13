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
        'CSN2':"#0000FF",\
        'CSN3':"#FF0000",\
        'CSN4':"#6495ED",\
        'CSN5':"#56341F",\
        'CSN6':"#40E0D0",\
        'CSN7':"#EECC00",\
        'CSN8':"#00FF00"
        }

totvol = 0.0
for p in prots:
    totvol += vol[p]

runCommand('set bgcolor white')
i=0

#runCommand('open ../CSN.pdb')
#i += 1

# Whole LPD
runCommand('open LPD_Whole.mrc')
runCommand('color #000000 #'+str(i))
runCommand('volume #'+str(i)+' step 1 transparency 0.50 encloseVolume ' + str(totvol))
i+=1    


# Read components
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1')
    runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' style bold size 52 xpos .1 ypos ' + str( 0.7 - i / 20.0)) 
    i += 1
    
# Read Sample A, then B
for p in prots:
    runCommand('open Sample_A/LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1')
    runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' hide')
    i += 1

for p in prots:
    runCommand('open Sample_B/LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' encloseVolume '+str(vol[p]/0.7)+' step 1')
    runCommand('volume #'+str(i)+' transparency 0.5')
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' hide') 
    i += 1

runCommand('open CSN.pdb')
runCommand('del #'+str(i)+':.I-Z')

for c in prots:
    runCommand('color '+col[c]+' #'+str(i)+':.'+chains[c])

runCommand('del #'+str(i)+':.G')
i += 1
runCommand('open CSN7B.pdb')
runCommand('color #EECC00 #'+str(i))
runCommand('center')

runCommand('windowsize 1024 1024')
runCommand('turn x 1.0 90')
runCommand('wait 10')
runCommand('move x 64')
runCommand('wait 10')
runCommand('scale 2.5')
